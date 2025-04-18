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

def convert_absorption_to_json(df: pd.DataFrame) -> dict:
    row = df.iloc[0].to_dict()

    return {
        "Absorption": {
            "description": "The absorption characteristics of a drug, including permeability and bioavailability.",
            "data": {
                "caco2": {
                    "value": row.get("caco2"),
                    "description": (
                        "Caco-2 Permeability (log cm/s). > -5.15: excellent (green); otherwise: poor (red)."
                    )
                },
                "PAMPA": {
                    "value": row.get("PAMPA"),
                    "description": (
                        "PAMPA logPeff. 0–0.3: excellent (green); 0.3–0.7: medium (yellow); 0.7–1.0: poor (red)."
                    )
                },
                "MDCK": {
                    "value": row.get("MDCK"),
                    "description": (
                        "MDCK Permeability (cm/s). >2 × 10⁻⁶ cm/s: excellent (green); otherwise: poor (red)."
                    )
                },
                "pgp_inh": {
                    "value": row.get("pgp_inh"),
                    "description": (
                        "P-glycoprotein inhibitor probability. 0–0.3: excellent; 0.3–0.7: medium; 0.7–1.0: poor."
                    )
                },
                "pgp_sub": {
                    "value": row.get("pgp_sub"),
                    "description": (
                        "P-glycoprotein substrate probability. 0–0.3: excellent; 0.3–0.7: medium; 0.7–1.0: poor."
                    )
                },
                "hia": {
                    "value": row.get("hia"),
                    "description": (
                        "Human intestinal absorption (HIA+) probability. 0–0.3: excellent; 0.3–0.7: medium; 0.7–1.0: poor."
                    )
                },
                "f20": {
                    "value": row.get("f20"),
                    "description": (
                        "Bioavailability < 20% probability (F20%+). 0–0.3: excellent; 0.3–0.7: medium; 0.7–1.0: poor."
                    )
                },
                "f30": {
                    "value": row.get("f30"),
                    "description": (
                        "Bioavailability < 30% probability (F30%+). 0–0.3: excellent; 0.3–0.7: medium; 0.7–1.0: poor."
                    )
                },
                "f50": {
                    "value": row.get("f50"),
                    "description": (
                        "Bioavailability < 50% probability (F50%+). 0–0.3: excellent; 0.3–0.7: medium; 0.7–1.0: poor."
                    )
                }
            }
        }
    }

import random

def get_absorption_info(smiles_str: str) -> str:
    """模拟获取吸收性质的 JSON 数据，返回带描述的随机值结构。"""
    def rand_prob():
        return round(random.uniform(0.0, 1.0), 2)

    def rand_permeability():
        return round(random.uniform(-6.0, -4.5), 2)  # 比如 caco2

    def rand_pampa():
        return round(random.uniform(0.0, 1.0), 2)

    def rand_mdck():
        return round(random.uniform(1e-7, 5e-6), 8)

    json_data = {
        "Absorption": {
            "description": "The absorption characteristics of a drug, including permeability and bioavailability.",
            "data": {
                "caco2": {
                    "value": rand_permeability(),
                    "description": "Caco-2 Permeability (log cm/s). > -5.15: excellent (green); otherwise: poor (red)."
                },
                "PAMPA": {
                    "value": rand_pampa(),
                    "description": "PAMPA logPeff. 0–0.3: excellent (green); 0.3–0.7: medium (yellow); 0.7–1.0: poor (red)."
                },
                "MDCK": {
                    "value": rand_mdck(),
                    "description": "MDCK Permeability (cm/s). >2 × 10⁻⁶ cm/s: excellent (green); otherwise: poor (red)."
                },
                "pgp_inh": {
                    "value": rand_prob(),
                    "description": "P-glycoprotein inhibitor probability. 0–0.3: excellent; 0.3–0.7: medium; 0.7–1.0: poor."
                },
                "pgp_sub": {
                    "value": rand_prob(),
                    "description": "P-glycoprotein substrate probability. 0–0.3: excellent; 0.3–0.7: medium; 0.7–1.0: poor."
                },
                "hia": {
                    "value": rand_prob(),
                    "description": "Human intestinal absorption (HIA+) probability. 0–0.3: excellent; 0.3–0.7: medium; 0.7–1.0: poor."
                },
                "f20": {
                    "value": rand_prob(),
                    "description": "Bioavailability < 20% probability (F20%+). 0–0.3: excellent; 0.3–0.7: medium; 0.7–1.0: poor."
                },
                "f30": {
                    "value": rand_prob(),
                    "description": "Bioavailability < 30% probability (F30%+). 0–0.3: excellent; 0.3–0.7: medium; 0.7–1.0: poor."
                },
                "f50": {
                    "value": rand_prob(),
                    "description": "Bioavailability < 50% probability (F50%+). 0–0.3: excellent; 0.3–0.7: medium; 0.7–1.0: poor."
                }
            }
        }
    }
    return json.dumps(json_data, ensure_ascii=False)


# def get_absorption_info(smiles_str: str) -> str:
    """查询药物的吸收性质（ADMET 中的 A，Absorption），包括肠道吸收、血脑屏障穿透、P-gp 作用等，评估药物从给药部位进入血液循环过程中的行为。提交 SMILES 表达式，返回 JSON 格式输出。"""
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
        json_data = convert_absorption_to_json(df)
        return json.dumps(json_data, ensure_ascii=False)
        

    except Exception as e:
        print("❌ 错误：", e)
        return {"error": str(e)}

    finally:
        driver.quit()

@tool
def absorption_agent() -> str:
    """查询环境中药物的吸收性质（ADMET 中的 A，Absorption），包括肠道吸收、血脑屏障穿透、P-gp 作用等，评估药物从给药部位进入血液循环过程中的行为。针对session里面的molecule，返回一个JSON，以SMILES为key。"""
    molecule_list = st.session_state.get("molecule", [])
    results = {}

    for molecule in molecule_list:
        smiles_str = Chem.MolToSmiles(molecule)
        for attempt in range(2):  # 最多尝试 2 次
            try:
                result_json = json.loads(get_absorption_info(smiles_str))
                results[smiles_str] = result_json
                break  # 成功就跳出 retry 循环
            except Exception as e:
                if attempt == 1:
                    # 第二次也失败了，才真正记录错误
                    results[smiles_str] = {
                        "Absorption": {
                            "error": f"查询失败（已重试）: {str(e)}"
                        }
                    }

    return json.dumps(results, indent=2, ensure_ascii=False)


if __name__ == "__main__":
    smiles = "CCO"
    result = absorption_agent(smiles)
    print(result)