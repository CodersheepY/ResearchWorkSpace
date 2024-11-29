from pymongo import MongoClient

mongodb_uri = "mongodb://localhost:27017"

try:
    client = MongoClient(mongodb_uri)
    print("Connection successful!")
    print(client.server_info())
except Exception as e:
    print(f"Connection failed: {e}")
