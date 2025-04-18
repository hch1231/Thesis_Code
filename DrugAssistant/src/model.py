from rdkit import Chem


def mock_model(smile, type):
    try:
        molecule = Chem.MolFromSmiles(smile)
        return [molecule]
    except Exception as e:
        return f"创建过程遇到错误：{e}"

if __name__ == "__main__":
    mock_model("CCCCN1CCCCC1C(=O)NC1=C(C)C=CC=C1C",84)