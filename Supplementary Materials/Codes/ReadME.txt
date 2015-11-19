Following are the steps of applying our complex-detection framework: (Kindly note that the codes are written in Python)

1- Use the BicAT tool to bicluster the gene expression dataset, using various biclustering algorithms.
2- Run the "Extract and Filter_Cluster_PPIs" code to extract and filter the PPIs corresponding to each bicluster. Please specify the file directories, the initial PPI file and the exported biclustering file from the BicAT tool.
3- Apply PE filtering to clean the generated biclusters PPIs using the "Apply PE_EMH_CODE" python code. Please specify the corresponding file directories. The codes of the PE filtering method, PE_EMH.py and PEWCC_EMH.py should be in the same directory.
4- Apply the complex-detection algorithms on each bicluster PPIs.
5- Merge and Filter the generated sets of detected complexes using the "Merge Detected Complexes" and Filter Detected Complexes" codes consecutively. Please specify the file directories.
6- Refine the total set of detected complexes using results using the "Final Merging Detected Complexes" code. Please specify the file directories.
7- Evaluate the results against the reference catalogue, here CYC2008, using the evaluation codes.
