import streamlit as st
from langchain import hub
from dotenv import load_dotenv
from langchain.agents import AgentExecutor, create_tool_calling_agent
from langchain_core.runnables.history import RunnableWithMessageHistory
from langchain_community.chat_message_histories import StreamlitChatMessageHistory
from langchain_community.callbacks.streamlit import StreamlitCallbackHandler
import os
import logging
import sys
import json
import re

from src.agent_tools import *
from src.call_llm import client
from langchain.prompts import (
    ChatPromptTemplate,
    PromptTemplate,
)
from langchain.chains import LLMChain

# ============================ 环境与日志配置 ============================
# load_dotenv("/home/hch/DrugAgent/.env")

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s - %(levelname)s - %(message)s",
    handlers=[
        # 保存所有命令行输出
        logging.FileHandler("agent_outputs.log", encoding="utf-8"),
        # 同时打印到终端
        logging.StreamHandler(sys.stdout),
    ],
)
logger = logging.getLogger(__name__)
# ======================================================================


def get_chat_history_text(history):
    """将 StreamlitChatMessageHistory 转换为纯文本。"""
    messages = history.messages
    history_text = ""
    for msg in messages:
        if msg.type == "human":
            history_text += f"用户：{msg.content}\n"
        elif msg.type == "ai":
            history_text += f"助手：{msg.content}\n"
    return history_text


@st.cache_resource
def load_prompt():
    return hub.pull("hwchase17/structured-chat-agent")


def single_mol():
    """Streamlit 入口函数。"""

    # -------------------- 工具 & Prompt --------------------
    tools = [
        create_mol_from_smile,
        visualization,
        physicochem_agent,
        medicinal_chemistry_agent,
        toxicology_agent,
        metabolism_agent,
        absorption_agent,
        distribution_agent,
        excretion_agent,
    ]

    agent_prompt = load_prompt()

    # 协调 Agent 的提示词
    coordinator_prompt = ChatPromptTemplate.from_template(
        open("/prompt/coordinator_prompt.txt").read()
    )
    coordinator_agent = LLMChain(llm=client, prompt=coordinator_prompt)

    # 任务分解 Agent 的提示词
    decomposer_prompt = ChatPromptTemplate.from_template(
        open("/prompt/decomposer_prompt.txt").read()
    )
    decomposer_agent = LLMChain(llm=client, prompt=decomposer_prompt)

    # 综合 Agent 的提示词
    synthesizer_prompt = ChatPromptTemplate.from_template(
        open("/prompt/synthesizer_prompt.txt").read()
    )
    # -------------------- AgentExecutor 缓存 --------------------
    if "agent_executor" not in st.session_state:
        structured_agent = create_tool_calling_agent(
            llm=client, tools=tools, prompt=agent_prompt
        )
        st.session_state.agent_executor = AgentExecutor(
            agent=structured_agent,
            tools=tools,
            verbose=True,
            max_iterations=10,
            handle_parsing_errors=True,
            return_intermediate_steps=True,
        )
        logger.info("Structured Chat Agent 初始化完成。")

    # -------------------- 历史记录缓存 --------------------
    if "internal_history" not in st.session_state:
        st.session_state.internal_history = StreamlitChatMessageHistory(
            key="internal_history"
        )

    if "display_history" not in st.session_state:
        st.session_state.display_history = StreamlitChatMessageHistory(
            key="display_history"
        )

    # -------------------- RunnableWithMessageHistory --------------------
    agent_with_chat_history = RunnableWithMessageHistory(
        st.session_state.agent_executor,
        lambda session_id: st.session_state.internal_history,
        input_messages_key="input",
        history_messages_key="internal_history",
        internal_history=st.session_state.internal_history,
    )

    # -------------------- 回显历史消息 --------------------
    for msg in st.session_state.display_history.messages:
        st.chat_message(msg.type).write(msg.content)

    # -------------------- 主交互循环 --------------------
    if prompt := st.chat_input():
        # 用户输入
        st.chat_message("user").write(prompt)
        st_callback = StreamlitCallbackHandler(st.container())
        config = {"configurable": {"session_id": "any"}, "callbacks": [st_callback]}

        # ➤ Step 1: 协调智能体提取任务结构
        chat_history_text = get_chat_history_text(st.session_state.display_history)
        coord_result = coordinator_agent.invoke(
            {"request": prompt, "chat_history": chat_history_text}
        )
        logger.info("Coordinator Agent Output: %s", coord_result)

        # 将原始输出存储
        raw_output = coord_result["text"]
        json_text = re.search(r"\{.*\}", raw_output, re.DOTALL).group()
        data = json.loads(json_text)

        goal = data["goal"]
        is_complex = data["is_complex"]
        execution_notes = data["execution_notes"]

        st.session_state.display_history.add_user_message(prompt)

        # ➤ Step 2: 根据复杂度选择后续执行路径
        if is_complex:
            decomposer_output = decomposer_agent.invoke({"request": prompt})
            logger.info("Decomposer Agent Output: %s", decomposer_output)
            agent_input = decomposer_output["text"]
        else:
            agent_input = prompt

        # ➤ Step 3: 统一传给专家智能体执行
        response = agent_with_chat_history.invoke({"input": agent_input}, config)
        logger.info("Structured Agent Output: %s", response)

        # ➤ Step 4: 综合智能体分析
        refine_prompt = PromptTemplate.from_template(synthesizer_prompt)
        refine_agent = LLMChain(llm=client, prompt=refine_prompt)

        raw_summary = refine_agent.run({
            "request": prompt,
            "subtasks": json.dumps(decomposer_output, ensure_ascii=False),
            "expert_results": json.dumps(response.get("output"), ensure_ascii=False),
        })

        # 解析并分支处理
        try:
            result = json.loads(raw_summary)
        except json.JSONDecodeError:
            st.error("无法解析综合智能体的输出，请检查提示词或 LLM 响应。")
            result = {"is_solved": False, "unsolved_subtask": "解析错误", "final_summary": raw_summary}

        if not result.get("is_solved", False):
            unsolved = result.get("unsolved_subtask", "未知子任务")
            # 反馈给协调智能体
            coord_feedback = coordinator_agent.invoke({
                "request": unsolved,
                "chat_history": get_chat_history_text(st.session_state.display_history)
            })
            logger.info("Coordinator Feedback Output: %s", coord_feedback)

            st.chat_message("assistant").write(
                f"检测到未完成子任务：“{unsolved}”，已反馈给协调智能体，继续处理。"
            )
            st.warning(f"当前进展：{result.get('final_summary')}")
        else:
            final_text = result.get("final_summary", "")
            st.chat_message("assistant").write(final_text)
            st.session_state.display_history.add_ai_message(final_text)



if __name__ == "__main__":
    single_mol()
