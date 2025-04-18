from selenium import webdriver
from selenium.webdriver.firefox.options import Options
from selenium.webdriver.firefox.service import Service
from selenium.webdriver.common.by import By
from selenium.webdriver.support.ui import WebDriverWait
from selenium.webdriver.support import expected_conditions as EC
import os, time
import pandas as pd
from langchain_core.tools import tool
import json
import pandas as pd
import streamlit as st
from rdkit import Chem

def convert_toxicology_to_json(df: pd.DataFrame) -> dict:
    row = df.iloc[0].to_dict()

    def toxicity_field(name, key, description):
        return {
            "value": row.get(key),
            "description": description
        }

    return {
        "Toxicology": {
            "description": (
                "毒理学（Toxicology）模块评估候选药物的安全性，涵盖神经、肝、肾、血液、遗传毒性等方面，"
                "帮助识别可能引发毒副作用的分子结构，为药物开发中的毒性筛选和风险评估提供支持。"
            ),
            "data": {
                "Neurotoxicity": toxicity_field("Neurotoxicity", "Neurotoxicity-DI", "药物诱导的神经毒性，概率值0-1。0–0.3: 优秀（绿色）；0.3–0.7: 中等（黄色）；0.7–1: 差（红色）"),
                "Ototoxicity": toxicity_field("Ototoxicity", "Ototoxicity", "耳毒性，损伤听觉系统。范围：0–0.3: 优秀；0.3–0.7: 中等；0.7–1: 差"),
                "Hematotoxicity": toxicity_field("Hematotoxicity", "Hematotoxicity", "血液毒性，可能损伤造血系统。范围参考上同"),
                "Nephrotoxicity": toxicity_field("Nephrotoxicity", "Nephrotoxicity-DI", "肾毒性，评估对肾脏的伤害可能性"),
                "Genotoxicity": toxicity_field("Genotoxicity", "Genotoxicity", "遗传毒性，对DNA造成损害"),
                "RPMI-8226 Cytotoxicity": toxicity_field("RPMI-8226 Cytotoxicity", "RPMI-8226", "多发性骨髓瘤细胞毒性"),
                "A549 Cytotoxicity": toxicity_field("A549 Cytotoxicity", "A549", "肺癌细胞毒性"),
                "HEK293 Cytotoxicity": toxicity_field("HEK293 Cytotoxicity", "HEK293", "胚胎肾细胞毒性"),
                "hERG Blockers (10μM)": toxicity_field("hERG Blockers (10μM)", "hERG-10um", "心脏毒性标志物，阻断hERG通道（10μM条件）"),
                "hERG Blockers": toxicity_field("hERG Blockers", "hERG", "心律失常风险评估"),
                "Human Hepatotoxicity": toxicity_field("Human Hepatotoxicity", "H-HT", "人体肝毒性预测"),
                "DILI": toxicity_field("DILI", "DILI", "药物性肝损伤风险"),
                "AMES Mutagenicity": toxicity_field("AMES Mutagenicity", "Ames", "致突变性（AMES试验）"),
                "Rat Oral Acute Toxicity": toxicity_field("Rat Oral Acute Toxicity", "LD50_oral", "大鼠口服急性毒性，概率值越大风险越高"),
                "FDAMDD": toxicity_field("FDAMDD", "FDAMDD", "FDA推荐最大日剂量"),
                "Skin Sensitization": toxicity_field("Skin Sensitization", "SkinSen", "皮肤致敏性评估"),
                "Carcinogenicity": toxicity_field("Carcinogenicity", "Carcinogenicity", "致癌性风险"),
                "Eye Irritation / Corrosion": toxicity_field("Eye Irritation / Corrosion", "EI", "眼部刺激/腐蚀性"),
                "Respiratory Toxicity": toxicity_field("Respiratory Toxicity", "Respiratory", "呼吸毒性"),
                "BCF": toxicity_field("BCF", "BCF", "生物富集系数，越高表示越易富集在生物体内"),
                "IGC50": toxicity_field("IGC50", "IGC50", "四膜虫生长抑制浓度（越高越毒）"),
                "LC50FM": toxicity_field("LC50FM", "LC50FM", "小鱼致死浓度（96h）"),
                "LC50DM": toxicity_field("LC50DM", "LC50DM", "水蚤致死浓度（48h）"),
                "NR-AR": toxicity_field("NR-AR", "NR-AR", "雄激素受体激动剂/拮抗剂活性"),
                "NR-AR-LBD": toxicity_field("NR-AR-LBD", "NR-AR-LBD", "雄激素受体配体结合域结合活性"),
                "NR-AhR": toxicity_field("NR-AhR", "NR-AhR", "芳烃受体激动剂活性"),
                "NR-Aromatase": toxicity_field("NR-Aromatase", "NR-Aromatase", "芳香化酶抑制作用"),
                "NR-ER": toxicity_field("NR-ER", "NR-ER", "雌激素受体激动剂活性"),
                "NR-ER-LBD": toxicity_field("NR-ER-LBD", "NR-ER-LBD", "雌激素受体配体结合域活性"),
                "NR-PPAR-gamma": toxicity_field("NR-PPAR-gamma", "NR-PPAR-gamma", "PPARγ 受体调节活性"),
                "SR-ARE": toxicity_field("SR-ARE", "SR-ARE", "氧化应激通路响应活性"),
                "SR-ATAD5": toxicity_field("SR-ATAD5", "SR-ATAD5", "DNA损伤响应（ATAD5相关）"),
                "SR-HSE": toxicity_field("SR-HSE", "SR-HSE", "热激应答相关毒性"),
                "SR-MMP": toxicity_field("SR-MMP", "SR-MMP", "线粒体膜电位影响"),
                "SR-p53": toxicity_field("SR-p53", "SR-p53", "p53介导的细胞损伤响应"),
                "Acute Toxicity Rule": toxicity_field("Acute Toxicity Rule", "EC", "结构报警片段与急性毒性相关"),
                "Genotoxic Carcinogenicity Rule": toxicity_field("Genotoxic Carcinogenicity Rule", "Genotoxic_Carcinogenicity_Mutagenicity", "结构报警与基因毒性相关"),
                "NonGenotoxic Carcinogenicity Rule": toxicity_field("NonGenotoxic Carcinogenicity Rule", "NonGenotoxic_Carcinogenicity", "非基因毒性致癌风险结构片段"),
                "Skin Sensitization Rule": toxicity_field("Skin Sensitization Rule", "Skin_Sensitization", "结构报警与皮肤过敏性相关"),
                "Aquatic Toxicity Rule": toxicity_field("Aquatic Toxicity Rule", "Acute_Aquatic_Toxicity", "结构报警与水体毒性相关"),
                "NonBiodegradable Rule": toxicity_field("NonBiodegradable Rule", "NonBiodegradable", "不可降解结构片段"),
                "SureChEMBL Rule": toxicity_field("SureChEMBL Rule", "SureChEMBL", "不利于药物开发的结构警示")
            }
        }
    }


import random
import json

@tool
def get_toxicology_info(smiles_str: str) -> str:
    """模拟生成毒理学（Toxicology）性质的 JSON 数据，结构完整、数据合理、支持测试调试。"""

    def rand_prob():
        return round(random.uniform(0, 1), 2)

    def rand_alert():
        return random.choice([0, 1])

    def rand_toxic_value():
        return round(random.uniform(0, 10), 2)

    def rand_concentration():
        return round(random.uniform(0.1, 500), 2)

    keys = {
        "Neurotoxicity-DI": rand_prob(),
        "Ototoxicity": rand_prob(),
        "Hematotoxicity": rand_prob(),
        "Nephrotoxicity-DI": rand_prob(),
        "Genotoxicity": rand_prob(),
        "RPMI-8226": rand_prob(),
        "A549": rand_prob(),
        "HEK293": rand_prob(),
        "hERG-10um": rand_prob(),
        "hERG": rand_prob(),
        "H-HT": rand_prob(),
        "DILI": rand_prob(),
        "Ames": rand_alert(),
        "LD50_oral": rand_toxic_value(),
        "FDAMDD": rand_alert(),
        "SkinSen": rand_prob(),
        "Carcinogenicity": rand_alert(),
        "EI": rand_prob(),
        "Respiratory": rand_prob(),
        "BCF": rand_toxic_value(),
        "IGC50": rand_concentration(),
        "LC50FM": rand_concentration(),
        "LC50DM": rand_concentration(),
        "NR-AR": rand_prob(),
        "NR-AR-LBD": rand_prob(),
        "NR-AhR": rand_prob(),
        "NR-Aromatase": rand_prob(),
        "NR-ER": rand_prob(),
        "NR-ER-LBD": rand_prob(),
        "NR-PPAR-gamma": rand_prob(),
        "SR-ARE": rand_prob(),
        "SR-ATAD5": rand_prob(),
        "SR-HSE": rand_prob(),
        "SR-MMP": rand_prob(),
        "SR-p53": rand_prob(),
        "EC": rand_alert(),
        "Genotoxic_Carcinogenicity_Mutagenicity": rand_alert(),
        "NonGenotoxic_Carcinogenicity": rand_alert(),
        "Skin_Sensitization": rand_alert(),
        "Acute_Aquatic_Toxicity": rand_alert(),
        "NonBiodegradable": rand_alert(),
        "SureChEMBL": rand_alert()
    }

    def toxicity_field(key, description):
        return {
            "value": keys[key],
            "description": description
        }

    json_data = {
        "Toxicology": {
            "description": (
                "毒理学（Toxicology）模块评估候选药物的安全性，涵盖神经、肝、肾、血液、遗传毒性等方面，"
                "帮助识别可能引发毒副作用的分子结构，为药物开发中的毒性筛选和风险评估提供支持。"
            ),
            "data": {
                "Neurotoxicity": toxicity_field("Neurotoxicity-DI", "药物诱导的神经毒性，概率值0-1。0–0.3: 优秀（绿色）；0.3–0.7: 中等（黄色）；0.7–1: 差（红色）"),
                "Ototoxicity": toxicity_field("Ototoxicity", "耳毒性，损伤听觉系统。范围：0–0.3: 优秀；0.3–0.7: 中等；0.7–1: 差"),
                "Hematotoxicity": toxicity_field("Hematotoxicity", "血液毒性，可能损伤造血系统。范围参考上同"),
                "Nephrotoxicity": toxicity_field("Nephrotoxicity-DI", "肾毒性，评估对肾脏的伤害可能性"),
                "Genotoxicity": toxicity_field("Genotoxicity", "遗传毒性，对DNA造成损害"),
                "RPMI-8226 Cytotoxicity": toxicity_field("RPMI-8226", "多发性骨髓瘤细胞毒性"),
                "A549 Cytotoxicity": toxicity_field("A549", "肺癌细胞毒性"),
                "HEK293 Cytotoxicity": toxicity_field("HEK293", "胚胎肾细胞毒性"),
                "hERG Blockers (10μM)": toxicity_field("hERG-10um", "心脏毒性标志物，阻断hERG通道（10μM条件）"),
                "hERG Blockers": toxicity_field("hERG", "心律失常风险评估"),
                "Human Hepatotoxicity": toxicity_field("H-HT", "人体肝毒性预测"),
                "DILI": toxicity_field("DILI", "药物性肝损伤风险"),
                "AMES Mutagenicity": toxicity_field("Ames", "致突变性（AMES试验）"),
                "Rat Oral Acute Toxicity": toxicity_field("LD50_oral", "大鼠口服急性毒性，概率值越大风险越高"),
                "FDAMDD": toxicity_field("FDAMDD", "FDA推荐最大日剂量"),
                "Skin Sensitization": toxicity_field("SkinSen", "皮肤致敏性评估"),
                "Carcinogenicity": toxicity_field("Carcinogenicity", "致癌性风险"),
                "Eye Irritation / Corrosion": toxicity_field("EI", "眼部刺激/腐蚀性"),
                "Respiratory Toxicity": toxicity_field("Respiratory", "呼吸毒性"),
                "BCF": toxicity_field("BCF", "生物富集系数，越高表示越易富集在生物体内"),
                "IGC50": toxicity_field("IGC50", "四膜虫生长抑制浓度（越高越毒）"),
                "LC50FM": toxicity_field("LC50FM", "小鱼致死浓度（96h）"),
                "LC50DM": toxicity_field("LC50DM", "水蚤致死浓度（48h）"),
                "NR-AR": toxicity_field("NR-AR", "雄激素受体激动剂/拮抗剂活性"),
                "NR-AR-LBD": toxicity_field("NR-AR-LBD", "雄激素受体配体结合域结合活性"),
                "NR-AhR": toxicity_field("NR-AhR", "芳烃受体激动剂活性"),
                "NR-Aromatase": toxicity_field("NR-Aromatase", "芳香化酶抑制作用"),
                "NR-ER": toxicity_field("NR-ER", "雌激素受体激动剂活性"),
                "NR-ER-LBD": toxicity_field("NR-ER-LBD", "雌激素受体配体结合域活性"),
                "NR-PPAR-gamma": toxicity_field("NR-PPAR-gamma", "PPARγ 受体调节活性"),
                "SR-ARE": toxicity_field("SR-ARE", "氧化应激通路响应活性"),
                "SR-ATAD5": toxicity_field("SR-ATAD5", "DNA损伤响应（ATAD5相关）"),
                "SR-HSE": toxicity_field("SR-HSE", "热激应答相关毒性"),
                "SR-MMP": toxicity_field("SR-MMP", "线粒体膜电位影响"),
                "SR-p53": toxicity_field("SR-p53", "p53介导的细胞损伤响应"),
                "Acute Toxicity Rule": toxicity_field("EC", "结构报警片段与急性毒性相关"),
                "Genotoxic Carcinogenicity Rule": toxicity_field("Genotoxic_Carcinogenicity_Mutagenicity", "结构报警与基因毒性相关"),
                "NonGenotoxic Carcinogenicity Rule": toxicity_field("NonGenotoxic_Carcinogenicity", "非基因毒性致癌风险结构片段"),
                "Skin Sensitization Rule": toxicity_field("Skin_Sensitization", "结构报警与皮肤过敏性相关"),
                "Aquatic Toxicity Rule": toxicity_field("Acute_Aquatic_Toxicity", "结构报警与水体毒性相关"),
                "NonBiodegradable Rule": toxicity_field("NonBiodegradable", "不可降解结构片段"),
                "SureChEMBL Rule": toxicity_field("SureChEMBL", "不利于药物开发的结构警示")
            }
        }
    }

    return json.dumps(json_data, indent=2, ensure_ascii=False)


# def get_toxicology_info(smiles_str: str) -> str:
    """查询分子的毒性信息（ADMET 中的 T，Toxicology），涵盖神经毒性、肝毒性、致癌性、遗传毒性等多个系统毒性预测指标，用于药物早期筛选阶段的毒理安全性评估。提交 SMILES 表达式，返回 JSON 格式输出。  """
    # 路径配置（你可根据需要外部传入）
    gecko_path = "/home/hch/DrugAgent/geckodriver"
    firefox_path = "/home/hch/DrugAgent/firefox/firefox"
    download_dir = "/home/hch/DrugAgent/downloads"

    # Firefox 设置
    options = Options()
    options.add_argument("--headless")
    options.binary_location = firefox_path
    options.set_preference("browser.download.folderList", 2)
    options.set_preference("browser.download.dir", download_dir)
    options.set_preference("browser.helperApps.neverAsk.saveToDisk", "text/csv")
    options.set_preference("pdfjs.disabled", True)

    # 启动 Firefox
    service = Service(executable_path=gecko_path)
    driver = webdriver.Firefox(service=service, options=options)

    try:
        # 打开提交页面
        driver.get("https://admetlab3.scbdd.com/server/evaluation")
        WebDriverWait(driver, 10).until(EC.presence_of_element_located((By.ID, "smiles")))

        # 输入 SMILES
        input_box = driver.find_element(By.ID, "smiles")
        input_box.clear()
        input_box.send_keys(smiles_str)

        # 提交
        driver.find_element(By.ID, "submit-btn-1").click()

        # 等待页面跳转
        WebDriverWait(driver, 30).until(
            lambda d: d.current_url != "https://admetlab3.scbdd.com/server/evaluation"
        )

        # 等待按钮加载
        WebDriverWait(driver, 20).until(EC.presence_of_element_located((By.ID, "download_csv")))

        # 清理旧文件
        for f in os.listdir(download_dir):
            if f.endswith(".csv"):
                os.remove(os.path.join(download_dir, f))

        # 点击下载
        csv_button = driver.find_element(By.ID, "download_csv")
        driver.execute_script("arguments[0].click();", csv_button)

        # 等待文件生成
        timeout = 30
        start_time = time.time()
        csv_path = None
        while time.time() - start_time < timeout:
            files = [f for f in os.listdir(download_dir) if f.endswith(".csv")]
            if files:
                csv_path = os.path.join(download_dir, files[0])
                break
            time.sleep(1)

        if not csv_path:
            raise RuntimeError("❌ CSV 文件未下载成功")

        # 读取 CSV 转成 JSON
        df = pd.read_csv(csv_path)
        json_data = convert_toxicology_to_json(df)
        return json.dumps(json_data, ensure_ascii=False)
        

    except Exception as e:
        print("❌ 错误：", e)
        return {"error": str(e)}

    finally:
        driver.quit()
    
@tool
def toxicology_agent() -> str:
    """查询环境中分子的毒性信息（ADMET 中的 T，Toxicology），涵盖神经毒性、肝毒性、致癌性、遗传毒性等多个系统毒性预测指标，用于药物早期筛选阶段的毒理安全性评估。返回 JSON 格式输出。  """
    molecule_list = st.session_state.get("molecule", [])
    results = {}

    for molecule in molecule_list:
        smiles_str = Chem.MolToSmiles(molecule)
        for attempt in range(2):  # 最多尝试两次
            try:
                result_json = json.loads(get_toxicology_info(smiles_str))
                results[smiles_str] = result_json
                break  # 成功就跳出循环
            except Exception as e:
                if attempt == 1:
                    results[smiles_str] = {
                        "Toxicology": {
                            "error": f"查询失败（已重试）: {str(e)}"
                        }
                    }

    return json.dumps(results, indent=2, ensure_ascii=False)

if __name__ == "__main__":
    smiles = "CCO"
    result = toxicology_agent(smiles)
    print(result)