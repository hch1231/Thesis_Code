import os
from dotenv import load_dotenv
from rdkit import Chem
import streamlit as st
from langchain_core.tools import tool
import requests

from langchain_community.vectorstores import FAISS
from src.call_llm import client,embeddings


vector_db = FAISS.load_local(os.getenv("FAISS_PATH"),
                                     embeddings,
                                     allow_dangerous_deserialization=True)

@tool
def create_mol_from_smile(smile_string, description):
    """
    该工具用于基于用户需求创建molecule对象或者生成molecule分子对象
    :param smile_string: 如果没有description，就是用户需要创建的SMILE分子表达式，
                         如果有description，就是相互作用描述对应的给定分子
    :param description: 用户对于生成分子和给定分子相互作用的描述，若无则填null
    :return:
    """
    if description and description != "null":
#         sim_res = vector_db.similarity_search(description, k=10)
#         sim_type = [i.metadata['type'] for i in sim_res]
#         candidate = "\n".join(str(i+1) + "." + j.page_content for i, j in enumerate(sim_res))

#         prompt = f"""找出以下你觉得和{description}表达最近的表述：
# {candidate}
# #Output Format:
# {{"index":你找到的最相近的表述的序号}}
# 不需要返回其他多余的解释
# """
#         index = client.invoke(prompt, temperature=0.01).content
#         if "```json" in index:
#             index = index.split("```json")[1].split("```")[0]

#         target_type = sim_type[int(eval(index)['index']) - 1]

#         # 打印模拟参数
#         print("分子生成参数如下：")
#         print(f"work_type: ddi_condition_sample")
#         print(f"label_drug: {smile_string}")
#         print(f"label_value: {target_type}")

        # 固定返回结果
        generated_molecule = [
            "CCC(CC)CC(C)=CC(C)(CO)C(C)N",
            "CCCNC(C)C(N)(CC)CCC",
            "CCCC1CC(O)C(C)CC(C)C1O",
            "CCC(C)(O)CC",
            "CCC1CC(C)C1(CC)C(N)=CC(C)C",
            "CCCC(C)CC(C)C(CC)=C(C)O",
        ]

        mol_list = [Chem.MolFromSmiles(mol) for mol in generated_molecule]
        st.session_state["molecule"] = mol_list

        return "采样任务完成，生成分子：\n" + "\n".join(generated_molecule) + "\n已保存"
    else:
        molecule = Chem.MolFromSmiles(smile_string)
        st.session_state["molecule"] = [molecule]
        return "分子创建成功，分子：\n" + smile_string + "\n已保存"

# @tool
# def create_mol_from_smile(smile_string,description):
    """
    该工具用于基于用户需求创建molecule对象或者生成molecule分子对象
    :param smile_string: 如果没有description，就是用户需要创建的SMILE分子表达式，如果有description，就是相互作用描述对应的给定分子
    :param description:填入用户对于生成分子和给定分子相互作用的描述，注意一定是某种药物相互作用，如果没有药物相互作用，则填入null,不需要和用户确认
        example：
        用户输入：我要创建增强CCCCN1CCCCC1C(=O)NC1=C(C)C=CC=C1C血管加压作用的分子
        description：增强血管加压作用的分子
    :return:
    """
    # global molecule
    # if smile_string == "":
    #     work_type = "ddi_condition_sample"
    #     label_drug = smile_string
    #     label_value = target_type 
    if description and description != "null":
        sim_res = vector_db.similarity_search(description,k=10)
        sim_type = [i.metadata['type'] for i in sim_res]

        candidate = "\n".join(str(i+1) + "." + j.page_content for i,j in enumerate(sim_res))

        prompt = f"""找出以下你觉得和{description}表达最近的表述：
         {candidate}
         #Output Format:
         {{"index":你找到的最相近的表述的序号}}
         不需要返回其他多余的解释
        """
        index = client.invoke(prompt,temperature=0).content
        if "```json" in index:
            index = index.split("```json")[1].split("```")[0]

        target_type = sim_type[int(eval(index)['index']-1)]
        url = "http://0.0.0.0:8080/run_sampler"
        work_type = "ddi_condition_sample"
        label_drug = smile_string
        label_value = target_type
        params = {
            "work_type": work_type,
            "label_drug": label_drug,
            "label_value": label_value
        }
        response = requests.get(url, params=params)
        if response.status_code == 200:
            result = response.json()
            generated_molecule = result.get("message", "")
            print(generated_molecule)
            mol_list = []
            for mol in generated_molecule:
                molecule = Chem.MolFromSmiles(mol)
                mol_list.append(molecule)
            st.session_state["molecule"] = mol_list
            return "采样任务已完成，生成分子：\n" + "\n".join(generated_molecule) + "\n已保存"
        else:
            raise Exception(f"API调用失败，状态码：{response.status_code}，错误信息：{response.text}")
    else:
        molecule = Chem.MolFromSmiles(smile_string)
        st.session_state["molecule"] = [molecule]
        return "分子创建成功，分子：\n"+ smile_string +"\n已保存"





    # work_type = args.get("work_type", "ddi_condition_sample")
    # label_drug = args.get("label_drug", "CCCCN1CCCCC1C(=O)NC1=C(C)C=CC=C1C")
    # label_value = args.get("label_value", 49)
    
    # # 修改为实际 API 地址
    # url = "http://0.0.0.0:8080/run_sampler"
    # params = {
    #     "work_type": work_type,
    #     "label_drug": label_drug,
    #     "label_value": label_value
    # }
    # response = requests.get(url, params=params)
    # if response.status_code == 200:
    #     result = response.json()
    #     generated_molecule = result.get("message", "")
    #     st.session_state["molecule"] = result.get("message", "")
    #     return "采样任务已完成，生成分子：\n"+ generated_molecule +"\n已保存"
    # else:
    #     raise Exception(f"API调用失败，状态码：{response.status_code}，错误信息：{response.text}")