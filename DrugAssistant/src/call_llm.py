
import os
from langchain_community.embeddings import OpenAIEmbeddings
from langchain_openai import ChatOpenAI
#导入环境变量
from dotenv import load_dotenv
load_dotenv()
print(os.getenv("OPENAI_KEY"))
print(os.getenv("BASE_URL"))

client = ChatOpenAI(openai_api_key="OPENAI_KEY",
                    base_url="BASE_URL",
                    temperature=0,
                    model_name="qwen2.5-14b-instruct",)

embeddings = OpenAIEmbeddings(
             openai_api_key='OPENAI_KEY',
        base_url="BASE_URL",
)

