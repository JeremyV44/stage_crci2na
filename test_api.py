import requests

def get_uniprot_data(gene_name, organism_id):
    base_url = "https://rest.uniprot.org/uniprotkb/search"
    query = f"gene_exact:{gene_name} AND organism_id:{organism_id}"
    params = {
        "query": query,
        "format": "json",
        "fields": "cc_function,keyword"    }
    response = requests.get(base_url, params=params)
    
    if response.status_code == 200:
        return response.json()
    else:
        print(f"Error: {response.status_code} - {response.text}")
        return None

def extract_keywords_info(uniprot_data):
    if uniprot_data.get('results'):
        entry = uniprot_data['results'][0]
        gene_info = {
            "biological_process": [],
            "cellular_component": [],
            "comment_value": ""
        }
        
        # Extraire les mots-clés
        for keyword in entry.get("keywords", []):
            if keyword.get("category") == "Biological process":
                gene_info["biological_process"].append(keyword.get("name", "N/A"))
            elif keyword.get("category") == "Cellular component":
                gene_info["cellular_component"].append(keyword.get("name", "N/A"))

        comments = entry.get("comments", [])
        for comment in comments:
            if comment.get("texts"):
                for text in comment["texts"]:
                    if "value" in text:
                        gene_info["comment_value"] = text["value"]
                        break  # Prendre seulement le premier texte
                break
        return gene_info
    
    return None


