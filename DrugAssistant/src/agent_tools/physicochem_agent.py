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

def convert_physicochem_to_json(df: pd.DataFrame) -> dict:
    
    row = df.iloc[0].to_dict()

    return {
        "Physicochemical Property": {
            "description": (
                "The physico-chemical properties of a system are thermodynamically defined by a set of "
                "macroscopic quantities accessible to experimental measurements and related to the laws "
                "of statistical mechanics governing its microscopic parts"
            ),
            "data": {
                "MW": {
                    "value": row.get("MW"),
                    "description": "Contain hydrogen atoms. Optimal:100~600, based on Drug-Like Soft rule."
                },
                "Vol": {
                    "value": row.get("Vol"),
                    "description": "Van der Waals volume."
                },
                "Dense": {
                    "value": row.get("Dense"),
                    "description": "Density = MW / Volume"
                },
                "nHA": {
                    "value": row.get("nHA"),
                    "description": "Number of hydrogen bond acceptors. Optimal:0~12"
                },
                "nHD": {
                    "value": row.get("nHD"),
                    "description": "Number of hydrogen bond donors. Optimal:0~7"
                },
                "nRot": {
                    "value": row.get("nRot"),
                    "description": "Number of rotatable bonds. Optimal:0~11"
                },
                "nRing": {
                    "value": row.get("nRing"),
                    "description": "Number of rings. Optimal:0~6"
                },
                "MaxRing": {
                    "value": row.get("MaxRing"),
                    "description": "Number of atoms in the biggest ring. Optimal:0~18"
                },
                "nHet": {
                    "value": row.get("nHet"),
                    "description": "Number of heteroatoms. Optimal:1~15"
                },
                "fChar": {
                    "value": row.get("fChar"),
                    "description": "Formal charge. Optimal:-4 ~ 4"
                },
                "nRig": {
                    "value": row.get("nRig"),
                    "description": "Number of rigid bonds. Optimal:0~30"
                },
                "Flex": {
                    "value": row.get("Flex"),
                    "description": "Flexibility = nRot / nRig"
                },
                "nStereo": {
                    "value": row.get("nStereo"),
                    "description": "Number of stereocenters. Optimal: ≤ 2"
                },
                "TPSA": {
                    "value": row.get("TPSA"),
                    "description": "Topological polar surface area. Optimal:0~140"
                },
                "logS": {
                    "value": row.get("logS"),
                    "description": "Aqueous solubility. Optimal: -4 ~ 0.5"
                },
                "logP": {
                    "value": row.get("logP"),
                    "description": "n-octanol/water distribution coefficient. Optimal: 0 ~ 3"
                },
                "logD": {
                    "value": row.get("logD"),
                    "description": "logD7.4, balance between lipophilicity and hydrophilicity"
                },
                "pka_acidic": {
                    "value": row.get("pka_acidic"),
                    "description": "Acid dissociation constant. Lower = stronger acid"
                },
                "pka_basic": {
                    "value": row.get("pka_basic"),
                    "description": "Base dissociation constant. Higher = stronger base"
                },
                "mp": {
                    "value": row.get("mp"),
                    "description": "Melting point in °C. <25: liquid, ≥25: solid"
                },
                "bp": {
                    "value": row.get("bp"),
                    "description": "Boiling point in °C. <25°C → gas"
                }
            }
        }
    }

import random
import json

def get_physicochem_info(smiles_str: str) -> str:
    """模拟生成物理化学性质（Physicochemical Property）的 JSON 数据，结构与真实数据一致。"""

    def rand(min_val, max_val, digits=2):
        return round(random.uniform(min_val, max_val), digits)

    json_data = {
        "Physicochemical Property": {
            "description": (
                "The physico-chemical properties of a system are thermodynamically defined by a set of "
                "macroscopic quantities accessible to experimental measurements and related to the laws "
                "of statistical mechanics governing its microscopic parts"
            ),
            "data": {
                "MW": {"value": rand(100, 600), "description": "Contain hydrogen atoms. Optimal:100~600, based on Drug-Like Soft rule."},
                "Vol": {"value": rand(150, 1000), "description": "Van der Waals volume."},
                "Dense": {"value": rand(0.5, 2.0), "description": "Density = MW / Volume"},
                "nHA": {"value": rand(0, 12, 0), "description": "Number of hydrogen bond acceptors. Optimal:0~12"},
                "nHD": {"value": rand(0, 7, 0), "description": "Number of hydrogen bond donors. Optimal:0~7"},
                "nRot": {"value": rand(0, 11, 0), "description": "Number of rotatable bonds. Optimal:0~11"},
                "nRing": {"value": rand(0, 6, 0), "description": "Number of rings. Optimal:0~6"},
                "MaxRing": {"value": rand(0, 18, 0), "description": "Number of atoms in the biggest ring. Optimal:0~18"},
                "nHet": {"value": rand(1, 15, 0), "description": "Number of heteroatoms. Optimal:1~15"},
                "fChar": {"value": rand(-4, 4, 0), "description": "Formal charge. Optimal:-4 ~ 4"},
                "nRig": {"value": rand(1, 30, 0), "description": "Number of rigid bonds. Optimal:0~30"},
                "Flex": {"value": rand(0.0, 3.0), "description": "Flexibility = nRot / nRig"},
                "nStereo": {"value": rand(0, 4, 0), "description": "Number of stereocenters. Optimal: ≤ 2"},
                "TPSA": {"value": rand(0, 140), "description": "Topological polar surface area. Optimal:0~140"},
                "logS": {"value": rand(-8, 1), "description": "Aqueous solubility. Optimal: -4 ~ 0.5"},
                "logP": {"value": rand(-2, 6), "description": "n-octanol/water distribution coefficient. Optimal: 0 ~ 3"},
                "logD": {"value": rand(-2, 5), "description": "logD7.4, balance between lipophilicity and hydrophilicity"},
                "pka_acidic": {"value": rand(1, 9), "description": "Acid dissociation constant. Lower = stronger acid"},
                "pka_basic": {"value": rand(5, 13), "description": "Base dissociation constant. Higher = stronger base"},
                "mp": {"value": rand(-50, 300), "description": "Melting point in °C. <25: liquid, ≥25: solid"},
                "bp": {"value": rand(30, 500), "description": "Boiling point in °C. <25°C → gas"}
            }
        }
    }

    return json.dumps(json_data, indent=2, ensure_ascii=False)


# def get_physicochem_info(smiles_str: str) -> str:
    """查询环境内药物的物理化学性质，一个系统的物理化学性质在热力学上由一组可通过实验测量的宏观量所定义，这些宏观量与支配其微观组成部分的统计力学规律密切相关。"""
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
        json_data = convert_physicochem_to_json(df)
        return json.dumps(json_data, ensure_ascii=False)
        

    except Exception as e:
        print("❌ 错误：", e)
        return {"error": str(e)}

    finally:
        driver.quit()

@tool
def physicochem_agent() -> str:
    """查询环境中药物的物理化学性质，一个系统的物理化学性质在热力学上由一组可通过实验测量的宏观量所定义，这些宏观量与支配其微观组成部分的统计力学规律密切相关。"""
    molecule_list = st.session_state.get("molecule", [])
    results = {}

    for molecule in molecule_list:
        smiles_str = Chem.MolToSmiles(molecule)
        for attempt in range(2):  # 最多尝试两次
            try:
                result_json = json.loads(get_physicochem_info(smiles_str))
                results[smiles_str] = result_json
                break
            except Exception as e:
                if attempt == 1:
                    results[smiles_str] = {
                        "Physicochemical Property": {
                            "error": f"查询失败（已重试）: {str(e)}"
                        }
                    }

    return json.dumps(results, indent=2, ensure_ascii=False)

if __name__ == "__main__":
    smiles = "CCO"
    result = get_physicochem_info(smiles)
    print(result)