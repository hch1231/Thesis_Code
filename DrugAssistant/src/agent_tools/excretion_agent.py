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

def convert_excretion_to_json(df: pd.DataFrame) -> dict:
    row = df.iloc[0].to_dict()

    return {
        "Excretion": {
            "description": (
                "排泄（Excretion）模块用于评估药物从体内清除的速度与机制，主要通过血浆清除率（CLplasma）"
                "和药物半衰期（T1/2）两个指标衡量。该信息对于确定给药频率、剂量设计和药代动力学建模至关重要。"
            ),
            "data": {
                "cl-plasma": {
                    "value": row.get("cl-plasma"),
                    "description": (
                        "血浆清除率（CLplasma），单位为 ml/min/kg，反映药物从血液中清除的速度。"
                        "范围解释：0–5 优秀（绿色）；5–15 中等（黄色）；>15 差（红色）。"
                    )
                },
                "t0.5": {
                    "value": row.get("t0.5"),
                    "description": (
                        "药物半衰期（T1/2），单位为小时，表示血药浓度降低一半所需时间。"
                        "范围解释：>8 小时为长半衰期（绿色）；1–8 小时为中等（黄色）；<1 小时为短半衰期（红色）。"
                    )
                }
            }
        }
    }

import random
import json

def get_excretion_info(smiles_str: str) -> str:
    """模拟获取药物的排泄性质，返回随机生成的 Excretion JSON 数据。"""

    def rand_clearance():  # 血浆清除率，单位 ml/min/kg
        return round(random.uniform(0, 30), 2)

    def rand_half_life():  # 半衰期，单位小时
        return round(random.uniform(0.2, 15), 2)

    json_data = {
        "Excretion": {
            "description": (
                "排泄（Excretion）模块用于评估药物从体内清除的速度与机制，主要通过血浆清除率（CLplasma）"
                "和药物半衰期（T1/2）两个指标衡量。该信息对于确定给药频率、剂量设计和药代动力学建模至关重要。"
            ),
            "data": {
                "cl-plasma": {
                    "value": rand_clearance(),
                    "description": (
                        "血浆清除率（CLplasma），单位为 ml/min/kg，反映药物从血液中清除的速度。"
                        "范围解释：0–5 优秀（绿色）；5–15 中等（黄色）；>15 差（红色）。"
                    )
                },
                "t0.5": {
                    "value": rand_half_life(),
                    "description": (
                        "药物半衰期（T1/2），单位为小时，表示血药浓度降低一半所需时间。"
                        "范围解释：>8 小时为长半衰期（绿色）；1–8 小时为中等（黄色）；<1 小时为短半衰期（红色）。"
                    )
                }
            }
        }
    }

    return json.dumps(json_data, indent=2, ensure_ascii=False)

# def get_excretion_info(smiles_str: str) -> str:
    """查询药物的排泄性质（ADMET 中的 E，Excretion），包括肾脏排泄、胆汁排泄、消化道排泄等指标，用于评估药物从体内清除的路径和速率，对药物代谢终末阶段的行为预测具有指导意义。提交 SMILES 表达式，返回 JSON 格式输出。  """
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
        json_data = convert_excretion_to_json(df)
        return json.dumps(json_data, ensure_ascii=False)
        

    except Exception as e:
        print("❌ 错误：", e)
        return {"error": str(e)}

    finally:
        driver.quit()

@tool
def excretion_agent() -> str:
    """查询环境中药物的排泄性质（ADMET 中的 E，Excretion），包括肾脏排泄、胆汁排泄、消化道排泄等指标，用于评估药物从体内清除的路径和速率，对药物代谢终末阶段的行为预测具有指导意义。返回 JSON 格式输出。  """
    molecule_list = st.session_state.get("molecule", [])
    results = {}

    for molecule in molecule_list:
        smiles_str = Chem.MolToSmiles(molecule)
        for attempt in range(2):  # 最多尝试两次
            try:
                result_json = json.loads(get_excretion_info(smiles_str))
                results[smiles_str] = result_json
                break  # 成功跳出循环
            except Exception as e:
                if attempt == 1:  # 第二次失败才报错
                    results[smiles_str] = {
                        "Excretion": {
                            "error": f"查询失败（已重试）: {str(e)}"
                        }
                    }

    return json.dumps(results, indent=2, ensure_ascii=False)


if __name__ == "__main__":
    smiles = "CCO"
    result = excretion_agent(smiles)
    print(result)