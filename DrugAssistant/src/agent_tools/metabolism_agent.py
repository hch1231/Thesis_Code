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

def convert_metabolism_to_json(df: pd.DataFrame) -> dict:
    row = df.iloc[0].to_dict()

    return {
        "Metabolism": {
            "description": (
                "代谢（Metabolism）模块评估候选药物在体内是否会成为多种细胞色素 P450 酶（CYPs）的底物或抑制剂，"
                "这些酶主要在肝脏中代谢药物。代谢特征决定了药物的稳定性、半衰期以及潜在的药物相互作用风险。"
                "此外，还包含人肝微粒体（HLM）稳定性预测，用于估计药物在体内的代谢清除速率。"
            ),
            "data": {
                "CYP1A2-inh": {
                    "value": row.get("CYP1A2-inh"),
                    "description": "预测为 CYP1A2 抑制剂的概率（0-1）。CYP1A2 为肝脏中代谢芳香族胺类化合物的重要酶，抑制可能导致药物相互作用。"
                },
                "CYP1A2-sub": {
                    "value": row.get("CYP1A2-sub"),
                    "description": "预测为 CYP1A2 底物的概率（0-1）。成为该酶底物意味着可能被其代谢，影响药物代谢速率。"
                },
                "CYP2C19-inh": {
                    "value": row.get("CYP2C19-inh"),
                    "description": "预测为 CYP2C19 抑制剂的概率。该酶对多种药物（如抗抑郁药、质子泵抑制剂）代谢重要，抑制可能影响血药浓度。"
                },
                "CYP2C19-sub": {
                    "value": row.get("CYP2C19-sub"),
                    "description": "预测为 CYP2C19 底物的概率。为底物表示该药可能受 CYP2C19 介导的代谢调控。"
                },
                "CYP2C9-inh": {
                    "value": row.get("CYP2C9-inh"),
                    "description": "预测为 CYP2C9 抑制剂的概率。CYP2C9 参与非甾体类抗炎药（NSAIDs）、华法林等药物的代谢。"
                },
                "CYP2C9-sub": {
                    "value": row.get("CYP2C9-sub"),
                    "description": "预测为 CYP2C9 底物的概率。为底物者需评估是否与其他药物发生代谢竞争。"
                },
                "CYP2D6-inh": {
                    "value": row.get("CYP2D6-inh"),
                    "description": "预测为 CYP2D6 抑制剂的概率。该酶具有基因多态性，对中枢神经类药物尤为重要。"
                },
                "CYP2D6-sub": {
                    "value": row.get("CYP2D6-sub"),
                    "description": "预测为 CYP2D6 底物的概率。可用于判断是否受代谢表型影响（快代谢/慢代谢者）。"
                },
                "CYP3A4-inh": {
                    "value": row.get("CYP3A4-inh"),
                    "description": "预测为 CYP3A4 抑制剂的概率。CYP3A4 是最主要的药物代谢酶，参与约 50% 药物代谢，抑制风险高。"
                },
                "CYP3A4-sub": {
                    "value": row.get("CYP3A4-sub"),
                    "description": "预测为 CYP3A4 底物的概率。为该酶底物者在与其他 CYP3A4 抑制剂联用时需关注药物暴露升高风险。"
                },
                "CYP2B6-inh": {
                    "value": row.get("CYP2B6-inh"),
                    "description": "预测为 CYP2B6 抑制剂的概率。该酶代谢一些中枢神经药物及麻醉剂。"
                },
                "CYP2B6-sub": {
                    "value": row.get("CYP2B6-sub"),
                    "description": "预测为 CYP2B6 底物的概率。对个体间差异敏感，可能影响药效或毒性。"
                },
                "CYP2C8-inh": {
                    "value": row.get("CYP2C8-inh"),
                    "description": "预测为 CYP2C8 抑制剂的概率。该酶参与特定抗癌药与糖尿病药物的代谢过程。"
                },
                "LM-human": {
                    "value": row.get("LM-human"),
                    "description": (
                        "人肝微粒体（HLM）稳定性预测值，代表药物代谢清除速率的概率（0-1）。"
                        "0–0.3 表示化合物稳定（half-life > 30 min），0.7–1 表示容易被肝脏快速代谢。"
                    )
                }
            }
        }
    }

import random
import json

def get_metabolism_info(smiles_str: str) -> str:
    """模拟生成代谢（Metabolism）性质的 JSON 数据，返回结构完整的随机值。"""

    def rand_prob():
        return round(random.uniform(0.0, 1.0), 2)

    json_data = {
        "Metabolism": {
            "description": (
                "代谢（Metabolism）模块评估候选药物在体内是否会成为多种细胞色素 P450 酶（CYPs）的底物或抑制剂，"
                "这些酶主要在肝脏中代谢药物。代谢特征决定了药物的稳定性、半衰期以及潜在的药物相互作用风险。"
                "此外，还包含人肝微粒体（HLM）稳定性预测，用于估计药物在体内的代谢清除速率。"
            ),
            "data": {
                "CYP1A2-inh": {"value": rand_prob(), "description": "预测为 CYP1A2 抑制剂的概率（0-1）。CYP1A2 为肝脏中代谢芳香族胺类化合物的重要酶，抑制可能导致药物相互作用。"},
                "CYP1A2-sub": {"value": rand_prob(), "description": "预测为 CYP1A2 底物的概率（0-1）。成为该酶底物意味着可能被其代谢，影响药物代谢速率。"},
                "CYP2C19-inh": {"value": rand_prob(), "description": "预测为 CYP2C19 抑制剂的概率。该酶对多种药物（如抗抑郁药、质子泵抑制剂）代谢重要，抑制可能影响血药浓度。"},
                "CYP2C19-sub": {"value": rand_prob(), "description": "预测为 CYP2C19 底物的概率。为底物表示该药可能受 CYP2C19 介导的代谢调控。"},
                "CYP2C9-inh": {"value": rand_prob(), "description": "预测为 CYP2C9 抑制剂的概率。CYP2C9 参与非甾体类抗炎药（NSAIDs）、华法林等药物的代谢。"},
                "CYP2C9-sub": {"value": rand_prob(), "description": "预测为 CYP2C9 底物的概率。为底物者需评估是否与其他药物发生代谢竞争。"},
                "CYP2D6-inh": {"value": rand_prob(), "description": "预测为 CYP2D6 抑制剂的概率。该酶具有基因多态性，对中枢神经类药物尤为重要。"},
                "CYP2D6-sub": {"value": rand_prob(), "description": "预测为 CYP2D6 底物的概率。可用于判断是否受代谢表型影响（快代谢/慢代谢者）。"},
                "CYP3A4-inh": {"value": rand_prob(), "description": "预测为 CYP3A4 抑制剂的概率。CYP3A4 是最主要的药物代谢酶，参与约 50% 药物代谢，抑制风险高。"},
                "CYP3A4-sub": {"value": rand_prob(), "description": "预测为 CYP3A4 底物的概率。为该酶底物者在与其他 CYP3A4 抑制剂联用时需关注药物暴露升高风险。"},
                "CYP2B6-inh": {"value": rand_prob(), "description": "预测为 CYP2B6 抑制剂的概率。该酶代谢一些中枢神经药物及麻醉剂。"},
                "CYP2B6-sub": {"value": rand_prob(), "description": "预测为 CYP2B6 底物的概率。对个体间差异敏感，可能影响药效或毒性。"},
                "CYP2C8-inh": {"value": rand_prob(), "description": "预测为 CYP2C8 抑制剂的概率。该酶参与特定抗癌药与糖尿病药物的代谢过程。"},
                "LM-human": {
                    "value": rand_prob(),
                    "description": (
                        "人肝微粒体（HLM）稳定性预测值，代表药物代谢清除速率的概率（0-1）。"
                        "0–0.3 表示化合物稳定（half-life > 30 min），0.7–1 表示容易被肝脏快速代谢。"
                    )
                }
            }
        }
    }

    return json.dumps(json_data, indent=2, ensure_ascii=False)


# def get_metabolism_info(smiles_str: str) -> str:
    """查询药物的代谢性质（ADMET 中的 M，Metabolism），包括酶代谢途径、CYP450 相关代谢、代谢稳定性等，用于评估药物在体内的分解与生物转化情况。提交 SMILES 表达式，返回 JSON 格式输出。"""
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
        json_data = convert_metabolism_to_json(df)
        return json.dumps(json_data, ensure_ascii=False)
        

    except Exception as e:
        print("❌ 错误：", e)
        return {"error": str(e)}

    finally:
        driver.quit()

@tool
def metabolism_agent() -> str:
    """查询环境中药物的代谢性质（ADMET 中的 M，Metabolism），包括酶代谢途径、CYP450 相关代谢、代谢稳定性等，用于评估药物在体内的分解与生物转化情况。返回 JSON 格式输出。"""
    molecule_list = st.session_state.get("molecule", [])
    results = {}

    for molecule in molecule_list:
        smiles_str = Chem.MolToSmiles(molecule)
        for attempt in range(2):  # 最多尝试两次
            try:
                result_json = json.loads(get_metabolism_info(smiles_str))
                results[smiles_str] = result_json
                break  # 成功则跳出循环
            except Exception as e:
                if attempt == 1:  # 第二次也失败了才报错
                    results[smiles_str] = {
                        "Metabolism": {
                            "error": f"查询失败（已重试）: {str(e)}"
                        }
                    }

    return json.dumps(results, indent=2, ensure_ascii=False)


if __name__ == "__main__":
    smiles = "CCO"
    result = metabolism_agent(smiles)
    print(result)