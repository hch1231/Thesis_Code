import requests

def predict_admet_single(smiles, verify_cert=True):
    """
    根据输入的 SMILES，先调用分子清洗接口，再调用单分子 ADMET 预测接口，
    返回预测结果数据。
    
    参数:
        smiles: 字符串，表示分子的 SMILES 表达式。
        verify_cert: 是否验证 HTTPS 证书，默认为 True，可根据需要设为 False。
        
    返回:
        ADMET 预测结果的 JSON 数据；如果清洗后没有有效分子，则返回错误提示。
    """
    base_url = "https://admetlab3.scbdd.com"
    
    # Step 1: 分子清洗调用 /api/washmol
    wash_url = f"{base_url}/api/washmol"
    wash_payload = {
        "SMILES": smiles,
        "feature": False,
        "uncertain": False
    }
    try:
        wash_response = requests.post(wash_url, json=wash_payload, verify=verify_cert)
        wash_response.raise_for_status()
    except Exception as e:
        raise Exception("分子清洗接口调用失败: " + str(e))
    
    wash_data = wash_response.json()
    print("分子清洗返回结果:", wash_data)
    
    # 根据返回数据格式判断清洗结果
    cleaned = None
    if isinstance(wash_data.get("data"), list):
        for s in wash_data["data"]:
            if "invalid" not in s.lower():
                cleaned = s
                break
    elif isinstance(wash_data.get("data"), str):
        if "invalid" not in wash_data.get("data").lower():
            cleaned = wash_data.get("data")
    else:
        raise Exception("分子清洗返回数据格式未知")
    
    if cleaned is None:
        return {"error": "分子清洗后没有有效分子"}
    
    # Step 2: 单分子 ADMET 预测调用 /api/single/admet
    admet_url = f"{base_url}/api/single/admet"
    admet_payload = {
        "SMILES": cleaned,
        "feature": True  # 如果需要计算更多描述符，可以将此参数设置为 True
    }
    try:
        admet_response = requests.post(admet_url, json=admet_payload, verify=verify_cert)
        admet_response.raise_for_status()
    except Exception as e:
        raise Exception("ADMET 预测接口调用失败: " + str(e))
    
    admet_result = admet_response.json()
    return admet_result

# 示例调用
if __name__ == '__main__':
    # 示例输入分子，确保输入的是一个有效的 SMILES 字符串
    test_smiles = "CC(C)OC(=O)CC(=O)CSC1=C(C=C2CCCC2=N1)C#N"
    try:
        result = predict_admet_single(test_smiles, verify_cert=False)
        print("预测结果:", result)
    except Exception as ex:
        print("错误信息:", ex)
