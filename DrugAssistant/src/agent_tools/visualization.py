import streamlit as st
from rdkit import Chem
from rdkit.Chem import Draw
from langchain_core.tools import tool


@tool
def visualization():
    """
    该工具用于可视化环境中的分子
    :return:
    """

    def st_show(molecule):
        img = Draw.MolToImage(molecule, size=(300, 300))
        smile = Chem.MolToSmiles(molecule)
        st.image(img, caption=f"{smile} Visualization")

    molecule = st.session_state["molecule"]

    if molecule:
        if not isinstance(molecule, list):
            st_show(molecule)

        else:
            for mol in molecule:
                st_show(mol)

        return "分子已成功可视化"
    else:
        st.error("无法解析 SMILES 表达式，生成分子失败。")