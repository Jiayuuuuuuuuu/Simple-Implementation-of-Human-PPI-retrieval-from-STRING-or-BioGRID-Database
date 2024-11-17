import streamlit as st
import requests
import pandas as pd
import networkx as nx
import matplotlib.pyplot as plt

#2(a)-Retrieve human PPI data from BioGRID
def retrieve_ppi_biogrid(target_protein):
    biogrid_url = "https://webservice.thebiogrid.org/interactions"
    params = {
        "accessKey": "a34c45e0df868a763cbde7aeb22d6331", 
        "format": "json",
        "searchNames": True,
        "geneList": target_protein,
        "organism": 9606,  # Human
        "searchbiogridids": True,
        "includeInteractors": True
    }
    
    response = requests.get(biogrid_url, params=params)
    network = response.json()
    
    network_df = pd.DataFrame.from_dict(network, orient='index')
    return network_df

#2(b)-Retrieve human PPI data from SPRING
def retrieve_ppi_string(target_protein):
    string_url = "https://string-db.org/api/json/network"
    params = {
        "identifiers": target_protein,
        "species": 9606  # Human
    }
    
    response = requests.get(string_url, params=params)
    
    network = response.json()
    network_df = pd.json_normalize(network)
    
    return network_df

#3-Create network graph
def generate_network(dataframe):
    global database
    if database == "STRING":
        column_a, column_b = "preferredName_A", "preferredName_B"
    else: 
        column_a, column_b = "OFFICIAL_SYMBOL_A", "OFFICIAL_SYMBOL_B"
    
    network_graph = nx.from_pandas_edgelist(dataframe, column_a, column_b)
    return network_graph

#4-Retrieve all the network centrality measures
def get_centralities(network_graph):
    centrality_measures = [
        ("Degree Centrality", nx.degree_centrality(network_graph)),
        ("Betweenness Centrality", nx.betweenness_centrality(network_graph)),
        ("Closeness Centrality", nx.closeness_centrality(network_graph)),
        ("Eigenvector Centrality", nx.eigenvector_centrality(network_graph, max_iter=1000)),
        ("PageRank", nx.pagerank(network_graph))
    ]
    return centrality_measures

st.title("Lab 2 - Ong Jia Yu")
st.header("Simple Implementation of Human PPI retrieval from STRING or BioGRID Database")
    
col1, col2 = st.columns(2)
with col1:
    protein_id = st.text_input("Enter Protein ID:")
with col2:
    database = st.selectbox("Select Database:", ["STRING", "BioGRID"])

if st.button("Retrieve"):
    if database == "STRING":
        df = retrieve_ppi_string(protein_id)
    else:
        df = retrieve_ppi_biogrid(protein_id)
    
    network_graph = generate_network(df)
        
    col1, col2 = st.columns(2)
    
    with col1:
        st.subheader("PPI Data Information")
        st.dataframe(df)
        st.subheader("Description of the network")
        st.write(f"Number of proteins(nodes) in network: {network_graph.number_of_nodes()}")
        st.write(f"Number of interactions(edges) in network: {network_graph.number_of_edges()}")
        
        centrality_measures = get_centralities(network_graph)
        degree_centrality = centrality_measures[0][1]
        highest_node = max(degree_centrality, key=degree_centrality.get)
        
        plt.figure(figsize=(10, 7))
        slayout = nx.spring_layout(network_graph, seed=123)
        nx.draw(network_graph, slayout, with_labels=True, node_size=1100, node_color='lightblue')
        nx.draw_networkx_nodes(network_graph, slayout, nodelist=[highest_node], node_size=1100, node_color='yellow')
        st.pyplot(plt)
    
    with col2:
        st.subheader("Centrality Measures")
        for measure_name, measure_values in centrality_measures:
            st.write(f"\n**{measure_name}**")

            df_measure = pd.DataFrame.from_dict(measure_values, 
                                              orient='index', 
                                              columns=['Centrality'])
            df_measure = df_measure.sort_values('Centrality', 
                                              ascending=False)
            st.dataframe(df_measure)