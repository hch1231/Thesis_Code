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

def convert_distribution_to_json(df: pd.DataFrame) -> dict:
    row = df.iloc[0].to_dict()

    return {
        "Distribution": {
            "description": "Distribution-related properties of the drug, including protein binding, permeability, and transporters.",
            "data": {
                "BCRP": {
                    "value": row.get("BCRP"),
                    "description": (
                        "Breast cancer resistance protein inhibitor probability. "
                        "0–0.3: poor (red); 0.3–0.7: medium (yellow); 0.7–1.0: excellent (green)."
                    )
                },
                "OATP1B1": {
                    "value": row.get("OATP1B1"),
                    "description": (
                        "OATP1B1 inhibitor probability. 0–0.3: excellent (green); 0.3–0.7: medium; 0.7–1.0: poor (red)."
                    )
                },
                "OATP1B3": {
                    "value": row.get("OATP1B3"),
                    "description": (
                        "OATP1B3 inhibitor probability. 0–0.3: excellent (green); 0.3–0.7: medium; 0.7–1.0: poor (red)."
                    )
                },
                "BSEP": {
                    "value": row.get("BSEP"),
                    "description": (
                        "BSEP inhibitor probability. 0–0.3: excellent (green); 0.3–0.7: medium; 0.7–1.0: poor (red)."
                    )
                },
                "MRP1": {
                    "value": row.get("MRP1"),
                    "description": (
                        "MRP1 inhibitor probability. 0–0.3: excellent (green); 0.3–0.7: medium; 0.7–1.0: poor (red)."
                    )
                },
                "PPB": {
                    "value": row.get("PPB"),
                    "description": (
                        "Plasma Protein Binding (%). ≤ 90%: excellent (green); > 90%: poor (red)."
                    )
                },
                "logVDss": {
                    "value": row.get("logVDss"),
                    "description": (
                        "Volume of distribution (L/kg). 0.04–20: excellent (green); otherwise: poor (red)."
                    )
                },
                "BBB": {
                    "value": row.get("BBB"),
                    "description": (
                        "Blood-brain barrier penetration probability. 0–0.3: excellent (green); "
                        "0.3–0.7: medium (yellow); 0.7–1.0: poor (red)."
                    )
                },
                "Fu": {
                    "value": row.get("Fu"),
                    "description": (
                        "Fraction unbound in plasma (%). ≥ 5%: excellent (green); < 5%: poor (red)."
                    )
                }
            }
        }
    }

import random

def get_distribution_info(smiles_str: str) -> str:
    """模拟获取分布性质（Distribution）的 JSON 数据，随机生成带描述的结构。"""

    def rand_prob():
        return round(random.uniform(0.0, 1.0), 2)

    def rand_ppb():
        return round(random.uniform(70, 100), 1)  # Plasma protein binding %

    def rand_vdss():
        return round(random.uniform(0.01, 25.0), 2)  # logVDss 近似为 L/kg 范围

    def rand_fu():
        return round(random.uniform(0.1, 10), 2)  # Fraction unbound %

    json_data = {
        "Distribution": {
            "description": "Distribution-related properties of the drug, including protein binding, permeability, and transporters.",
            "data": {
                "BCRP": {
                    "value": rand_prob(),
                    "description": "Breast cancer resistance protein inhibitor probability. 0–0.3: poor (red); 0.3–0.7: medium (yellow); 0.7–1.0: excellent (green)."
                },
                "OATP1B1": {
                    "value": rand_prob(),
                    "description": "OATP1B1 inhibitor probability. 0–0.3: excellent (green); 0.3–0.7: medium; 0.7–1.0: poor (red)."
                },
                "OATP1B3": {
                    "value": rand_prob(),
                    "description": "OATP1B3 inhibitor probability. 0–0.3: excellent (green); 0.3–0.7: medium; 0.7–1.0: poor (red)."
                },
                "BSEP": {
                    "value": rand_prob(),
                    "description": "BSEP inhibitor probability. 0–0.3: excellent (green); 0.3–0.7: medium; 0.7–1.0: poor (red)."
                },
                "MRP1": {
                    "value": rand_prob(),
                    "description": "MRP1 inhibitor probability. 0–0.3: excellent (green); 0.3–0.7: medium; 0.7–1.0: poor (red)."
                },
                "PPB": {
                    "value": rand_ppb(),
                    "description": "Plasma Protein Binding (%). ≤ 90%: excellent (green); > 90%: poor (red)."
                },
                "logVDss": {
                    "value": rand_vdss(),
                    "description": "Volume of distribution (L/kg). 0.04–20: excellent (green); otherwise: poor (red)."
                },
                "BBB": {
                    "value": rand_prob(),
                    "description": "Blood-brain barrier penetration probability. 0–0.3: excellent (green); 0.3–0.7: medium (yellow); 0.7–1.0: poor (red)."
                },
                "Fu": {
                    "value": rand_fu(),
                    "description": "Fraction unbound in plasma (%). ≥ 5%: excellent (green); < 5%: poor (red)."
                }
            }
        }
    }

    return json.dumps(json_data, indent=2, ensure_ascii=False)


# def get_distribution_info(smiles_str: str) -> str:
    
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
        json_data = convert_distribution_to_json(df)
        return json.dumps(json_data, ensure_ascii=False)
        

    except Exception as e:
        print("❌ 错误：", e)
        return {"error": str(e)}

    finally:
        driver.quit()

@tool
def distribution_agent() -> str:
    """查询环境中药物的分布性质（ADMET 中的 D，Distribution），包括血浆蛋白结合率、组织分布、分布体积等，用于评估药物在体内不同组织间的转运与分布情况。返回 JSON 格式输出。  """
    molecule_list = st.session_state.get("molecule", [])
    results = {}

    for molecule in molecule_list:
        smiles_str = Chem.MolToSmiles(molecule)
        for attempt in range(2):  # 最多尝试两次
            try:
                result_json = json.loads(get_distribution_info(smiles_str))
                results[smiles_str] = result_json
                break
            except Exception as e:
                if attempt == 1:
                    results[smiles_str] = {
                        "Distribution": {
                            "error": f"查询失败（已重试）: {str(e)}"
                        }
                    }

    return json.dumps(results, indent=2, ensure_ascii=False)


if __name__ == "__main__":
    smiles = "CCO"
    result = distribution_agent(smiles)
    print(result)