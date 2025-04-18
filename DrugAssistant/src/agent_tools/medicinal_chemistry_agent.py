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
import streamlit as st
from rdkit import Chem

def convert_medicinal_chemistry_to_json(df: pd.DataFrame) -> dict:
    row = df.iloc[0].to_dict()

    return {
        "Medicinal Chemistry": {
            "description": "Drug-likeness and synthetic feasibility evaluations, including classic rules and expert-derived filters.",
            "data": {
                "QED": {
                    "value": row.get("QED"),
                    "description": (
                        "A measure of drug-likeness based on the concept of desirability. > 0.67: excellent (green); ≤ 0.67: poor (red)."
                    )
                },
                "Synth": {
                    "value": row.get("Synth"),
                    "description": (
                        "Synthetic accessibility score. ≤ 6: easy to synthesize (green); > 6: difficult (red)."
                    )
                },
                "gasa": {
                    "value": row.get("gasa"),
                    "description": (
                        "Graph-based synthetic accessibility score. 1: easy (green); 0: hard (red)."
                    )
                },
                "Fsp3": {
                    "value": row.get("Fsp3"),
                    "description": (
                        "Fraction of sp3 hybridized carbons. ≥ 0.42: good (green); < 0.42: poor (red)."
                    )
                },
                "MCE-18": {
                    "value": row.get("MCE-18"),
                    "description": (
                        "Medicinal chemistry evolution score. ≥ 45: high novelty (green); < 45: trivial (red)."
                    )
                },
                "Natural Product-likeness": {
                    "value": row.get("Natural Product-likeness"),
                    "description": (
                        "Natural product-likeness score (-5 ~ 5). Higher = more NP-like."
                    )
                },
                "Lipinski": {
                    "value": row.get("Lipinski"),
                    "description": (
                        "MW ≤ 500, logP ≤ 5, Hacc ≤ 10, Hdon ≤ 5. <2 violations: good (green); ≥2: poor (red)."
                    )
                },
                "Pfizer": {
                    "value": row.get("Pfizer"),
                    "description": (
                        "logP > 3, TPSA < 75. Both satisfied = poor (red); otherwise good."
                    )
                },
                "GSK": {
                    "value": row.get("GSK"),
                    "description": (
                        "MW ≤ 400, logP ≤ 4. 0 violations = good (green); otherwise poor."
                    )
                },
                "GoldenTriangle": {
                    "value": row.get("GoldenTriangle"),
                    "description": (
                        "200 ≤ MW ≤ 500; -2 ≤ logD ≤ 5. Fully met = excellent."
                    )
                },
                "PAINS": {
                    "value": row.get("PAINS"),
                    "description": (
                        "Pan-assay interference alerts (0 = safe). More alerts = more caution."
                    )
                },
                "Alarm_NMR": {
                    "value": row.get("Alarm_NMR"),
                    "description": (
                        "ALARM NMR: thiol-reactive substructures. 0 = good."
                    )
                },
                "BMS": {
                    "value": row.get("BMS"),
                    "description": (
                        "BMS rule: reactive substructures. 0 = preferred."
                    )
                },
                "Chelating": {
                    "value": row.get("Chelating"),
                    "description": (
                        "Chelating rule. 0 = non-chelator."
                    )
                },
                "Aggregators": {
                    "value": row.get("Aggregators"),
                    "description": (
                        "Colloidal aggregation probability. 0 = safe; >0 = caution."
                    )
                },
                "Fluc": {
                    "value": row.get("Fluc"),
                    "description": (
                        "FLuc inhibition. 0 = non-inhibitor; >0 = potential assay interference."
                    )
                },
                "Blue_fluorescence": {
                    "value": row.get("Blue_fluorescence"),
                    "description": (
                        "Blue fluorescence probability. 0 = clean."
                    )
                },
                "Green_fluorescence": {
                    "value": row.get("Green_fluorescence"),
                    "description": (
                        "Green fluorescence probability. 0 = clean."
                    )
                },
                "Reactive": {
                    "value": row.get("Reactive"),
                    "description": (
                        "Reactivity probability. 0 = stable."
                    )
                },
                "Promiscuous": {
                    "value": row.get("Promiscuous"),
                    "description": (
                        "Promiscuity probability. 0 = selective; >0 = potential off-targets."
                    )
                }
            }
        }
    }

import random
import json

def get_medicinal_chemistry_info(smiles_str: str) -> str:
    """模拟生成药物化学性质评估的 JSON 数据，返回结构完整、带描述的随机值。"""

    def rand_float(min_val, max_val, digits=2):
        return round(random.uniform(min_val, max_val), digits)

    def rand_int(min_val, max_val):
        return random.randint(min_val, max_val)

    def rand_binary():
        return random.choice([0, 1])

    json_data = {
        "Medicinal Chemistry": {
            "description": "Drug-likeness and synthetic feasibility evaluations, including classic rules and expert-derived filters.",
            "data": {
                "QED": {"value": rand_float(0.3, 0.95), "description": "A measure of drug-likeness based on the concept of desirability. > 0.67: excellent (green); ≤ 0.67: poor (red)."},
                "Synth": {"value": rand_float(2, 9), "description": "Synthetic accessibility score. ≤ 6: easy to synthesize (green); > 6: difficult (red)."},
                "gasa": {"value": rand_binary(), "description": "Graph-based synthetic accessibility score. 1: easy (green); 0: hard (red)."},
                "Fsp3": {"value": rand_float(0.0, 1.0), "description": "Fraction of sp3 hybridized carbons. ≥ 0.42: good (green); < 0.42: poor (red)."},
                "MCE-18": {"value": rand_int(10, 100), "description": "Medicinal chemistry evolution score. ≥ 45: high novelty (green); < 45: trivial (red)."},
                "Natural Product-likeness": {"value": rand_float(-5.0, 5.0), "description": "Natural product-likeness score (-5 ~ 5). Higher = more NP-like."},
                "Lipinski": {"value": rand_int(0, 4), "description": "MW ≤ 500, logP ≤ 5, Hacc ≤ 10, Hdon ≤ 5. <2 violations: good (green); ≥2: poor (red)."},
                "Pfizer": {"value": rand_binary(), "description": "logP > 3, TPSA < 75. Both satisfied = poor (red); otherwise good."},
                "GSK": {"value": rand_int(0, 2), "description": "MW ≤ 400, logP ≤ 4. 0 violations = good (green); otherwise poor."},
                "GoldenTriangle": {"value": rand_binary(), "description": "200 ≤ MW ≤ 500; -2 ≤ logD ≤ 5. Fully met = excellent."},
                "PAINS": {"value": rand_int(0, 2), "description": "Pan-assay interference alerts (0 = safe). More alerts = more caution."},
                "Alarm_NMR": {"value": rand_binary(), "description": "ALARM NMR: thiol-reactive substructures. 0 = good."},
                "BMS": {"value": rand_binary(), "description": "BMS rule: reactive substructures. 0 = preferred."},
                "Chelating": {"value": rand_binary(), "description": "Chelating rule. 0 = non-chelator."},
                "Aggregators": {"value": rand_binary(), "description": "Colloidal aggregation probability. 0 = safe; >0 = caution."},
                "Fluc": {"value": rand_binary(), "description": "FLuc inhibition. 0 = non-inhibitor; >0 = potential assay interference."},
                "Blue_fluorescence": {"value": rand_binary(), "description": "Blue fluorescence probability. 0 = clean."},
                "Green_fluorescence": {"value": rand_binary(), "description": "Green fluorescence probability. 0 = clean."},
                "Reactive": {"value": rand_binary(), "description": "Reactivity probability. 0 = stable."},
                "Promiscuous": {"value": rand_binary(), "description": "Promiscuity probability. 0 = selective; >0 = potential off-targets."}
            }
        }
    }

    return json.dumps(json_data, indent=2, ensure_ascii=False)


# def get_medicinal_chemistry_info(smiles_str: str) -> str:
    """查询化合物的成药性，药物化学涉及适用于治疗的新化学实体的识别、合成与开发。它还包括对现有药物、生物学性质及其定量构效关系（QSAR）的研究。"""
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

        json_data = convert_medicinal_chemistry_to_json(df)
        return json.dumps(json_data, ensure_ascii=False)
        

    except Exception as e:
        print("❌ 错误：", e)
        return {"error": str(e)}

    finally:
        driver.quit()

@tool
def medicinal_chemistry_agent() -> str:
    """查询环境中化合物的成药性，药物化学涉及适用于治疗的新化学实体的识别、合成与开发。它还包括对现有药物、生物学性质及其定量构效关系（QSAR）的研究。"""
    molecule_list = st.session_state.get("molecule", [])
    results = {}

    for molecule in molecule_list:
        smiles_str = Chem.MolToSmiles(molecule)
        for attempt in range(2):  # 最多尝试两次
            try:
                result_json = json.loads(get_medicinal_chemistry_info(smiles_str))
                results[smiles_str] = result_json
                break  # 成功就跳出重试
            except Exception as e:
                if attempt == 1:
                    results[smiles_str] = {
                        "Medicinal Chemistry": {
                            "error": f"查询失败（已重试）: {str(e)}"
                        }
                    }

    return json.dumps(results, indent=2, ensure_ascii=False)


if __name__ == "__main__":
    smiles = "CCO"
    result = get_medicinal_chemistry_info(smiles)
    print(result)