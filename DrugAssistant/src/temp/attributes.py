import requests
import time
import json
import sys
import pandas as pd
from rdkit import Chem
from rdkit.Chem import Descriptors
from langchain_core.tools import tool

def log(msg):
    print(msg)

def get_pubchem_properties(smiles):
    """查询 PubChem API 以获取分子性质"""
    url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/smiles/{smiles}/property/MolecularWeight,XLogP,HBondDonorCount,HBondAcceptorCount,TPSA/JSON"
    response = requests.get(url)
    
    if response.status_code == 200:
        data = response.json()
        props = data['PropertyTable']['Properties'][0]
        return {
            '分子量': float(props.get('MolecularWeight', 0)),
            'LogP': float(props.get('XLogP', 0)),
            '氢键供体': int(props.get('HBondDonorCount', 0)),
            '氢键受体': int(props.get('HBondAcceptorCount', 0)),
            'TPSA': float(props.get('TPSA', 0))
        }
    else:
        log("PubChem API 请求失败: " + str(response.status_code))
        return None

def calculate_properties(smiles):
    """计算额外的分子性质"""
    mol = Chem.MolFromSmiles(smiles)
    if mol:
        return {
            '可旋转键': int(Descriptors.NumRotatableBonds(mol)),
            '摩尔折射率': float(Descriptors.MolMR(mol))
        }
    return {'可旋转键': 0, '摩尔折射率': 0}

def evaluate_druggability(properties):
    """评估分子的成药性"""
    mw, logp, hbd, hba, tpsa = (properties['分子量'], properties['LogP'], 
                               properties['氢键供体'], properties['氢键受体'], properties['TPSA'])
    rb, mr = properties['可旋转键'], properties['摩尔折射率']
    
    lipinski = "符合 Lipinski 规则" if sum([mw > 500, logp > 5, hbd > 5, hba > 10, rb > 10]) < 2 else "可能口服利用度较差"
    veber = "符合 Veber 规则" if rb <= 10 and tpsa <= 140 else "可能吸收较差"
    ghose = "符合 Ghose 规则" if (160 <= mw <= 480 and 0.4 <= logp <= 5.6 and 2 <= hba <= 20 and 40 <= mr <= 130) else "可能不是药物样分子"
    egan = "符合 Egan 规则" if logp <= 5.88 and tpsa <= 131 else "可能吸收性较低"
    muegge = "符合 Muegge 规则" if (200 <= mw <= 600 and logp <= 5 and tpsa <= 150 and hba <= 10 and hbd <= 5 and rb <= 15) else "可能不符合广义药物相似性"
    
    return {
        'Lipinski 规则': lipinski,
        'Veber 规则': veber,
        'Ghose 规则': ghose,
        'Egan 规则': egan,
        'Muegge 规则': muegge
    }

def dict_to_markdown_table(data):
    md = "| Key | Value |\n| --- | --- |\n"
    for key, value in data.items():
        md += f"| {key} | {value} |\n"
    return md

@tool
def get_molecular_properties(smiles):
    """查询化合物的全部性质，包括 PubChem 和计算属性，并评估成药性"""
    log("查询 PubChem 数据...")
    pubchem_data = get_pubchem_properties(smiles)
    
    if not pubchem_data:
        return "无法从 PubChem 获取数据，请检查 SMILES 是否正确"
    
    log("计算额外的分子性质...")
    extra_properties = calculate_properties(smiles)
    all_properties = {**pubchem_data, **extra_properties}
    
    log("评估成药性...")
    druggability = evaluate_druggability(all_properties)
    
    result = {**all_properties, **druggability}
    result = dict_to_markdown_table(result)
    return result


def main():
    smiles_list = ["CN(C(=O)CN(CCO)CC(=O)N(C)C(C)(C)CC1=CC=CC=C1)C(C)(C)CC1=CC=CC=C1"]
    results = []
    
    for smiles in smiles_list:
        log(f"正在处理 SMILES: {smiles}")
        result = get_molecular_properties(smiles)
        if isinstance(result, dict):
            results.append(result)
        time.sleep(2)  # 避免过快请求 API
    
    df = pd.DataFrame(results)
    print(df.to_markdown(index=False))

if __name__ == "__main__":
    main()