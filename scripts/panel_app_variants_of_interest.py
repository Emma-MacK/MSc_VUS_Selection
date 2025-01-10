import json
import requests

with open("filtered_variants_genes_of_interest.csv") as file:
    entry_list = file.read()

split_list = list(entry_list.split("\n"))
#split_list = ["DDX3X", "SLC2A1"]
# print(split_list)
server = "https://panelapp.genomicsengland.co.uk/api/v1"

new_csv = ""

for entry in split_list:
    gene = entry.split(",")[1]
    # print(gene)
    entry_panel_full = entry.split(",")[11]
    entry_panel = entry_panel_full.split("-")[0]
    # print(entry_panel)
    gene = gene.upper()
    ext = "/genes/" + gene
    content_headers = {"Content-Type": "application/json"}
    r = requests.get(server+ext, headers=content_headers, timeout=120)
    result_panelapp = r.json()

    if int(result_panelapp["count"]) > 0:

        result_list=result_panelapp["results"]
        # print(result_list)
        # print("PanelApp entry detected for ", gene)
        txt_fill = ",NA"
        for i in range(0,len(result_list)):

            results_dict=result_list[i]
            panel=results_dict["panel"]
            # print(panel)
            if panel["name"] == entry_panel:

                if "Expert Review Green" in results_dict["evidence"]:
                    txt_fill = ",Green for " + panel["name"]

                elif "Expert Review Amber" in results_dict["evidence"]:
                    txt_fill = ",Amber for " + panel["name"]

                elif "Expert Review Red" in results_dict["evidence"]:
                    txt_fill = ",Red for " + panel["name"]

            elif entry_panel in panel["relevant_disorders"]:

                if "Expert Review Green" in results_dict["evidence"]:
                    txt_fill = "," + entry_panel + " is included in relevant disorders for panel " + panel["name"] + " and gene is Green for this panel"

                elif "Expert Review Amber" in results_dict["evidence"]:
                    txt_fill = "," + entry_panel + " is included in relevant disorders for panel " + panel["name"] + " and gene is Amber for this panel"

                elif "Expert Review Red" in results_dict["evidence"]:
                    txt_fill = "," + entry_panel + " is included in relevant disorders for panel " + panel["name"] + " and gene is Red for this panel"
        print(entry, txt_fill)

    else:
        print(entry,",panel_in Panelapp")

# get HGNC ids and print with gene list
