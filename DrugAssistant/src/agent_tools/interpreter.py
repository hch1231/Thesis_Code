import os
import yaml
import streamlit as st
from dotenv import load_dotenv
from langchain_core.tools import tool

_ = load_dotenv()
with open(os.getenv("PROMPT_PATH"), "r",encoding="utf-8") as file:
    prompt = yaml.load(file, Loader=yaml.FullLoader)

multi_mol_template = prompt['multi_mol_template']
# single_mol_template = prompt['single_mol_template']

@tool
def python_inter(py_code):
    """
    Default docstring for python_inter.
    """
    if "molecule" not in st.session_state:
        return "Error: No molecules loaded. Please upload a file first."

    global molecule
    molecule = st.session_state["molecule"]

    global get_physicochem_info
    get_physicochem_info = st.session_state["get_physicochem_info"]

    global get_medicinal_chemistry_info
    get_medicinal_chemistry_info = st.session_state["get_medicinal_chemistry_info"]

    global get_absorption_info
    get_absorption_info = st.session_state["get_absorption_info"]

    global get_distribution_info
    get_distribution_info = st.session_state["get_distribution_info"]
    
    global get_metabolism_info
    get_metabolism_info = st.session_state["get_metabolism_info"]

    global get_excretion_info
    get_excretion_info = st.session_state["get_excretion_info"]

    global get_toxicology_info
    get_toxicology_info = st.session_state["get_toxicology_info"]



    global_vars_before = set(globals().keys())
    try:
        exec(py_code, globals())
    except Exception as e:
        return f"代码执行时报错{e}"
    global_vars_after = set(globals().keys())
    new_vars = global_vars_after - global_vars_before
    # 若存在新变量
    if new_vars:
        result = {var: globals()[var] for var in new_vars}
        return str(result)
    # 若不存在新变量，即有可能是代码是表达式，也有可能代码对相同变量重复赋值
    else:
        try:
            # 尝试如果是表达式，则返回表达式运行结果
            return str(eval(py_code, globals()))
        # 若报错，则先测试是否是对相同变量重复赋值
        except Exception as e:
            try:
                exec(py_code, globals())
                return "已经顺利执行代码"
            except Exception as e:
                pass
            # 若不是重复赋值，则报错
            return f"代码执行时报错{e}"


DESCRIPTION_TEMPLATE = multi_mol_template
python_inter.description = DESCRIPTION_TEMPLATE
