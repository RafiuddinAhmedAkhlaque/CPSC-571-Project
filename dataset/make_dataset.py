import pandas as pd
import mysql.connector as sql
import csv


# Connect to Rfam database and get ncRNA sequence IDs we're interested in
rfam_db = sql.connect(
    host="mysql-rfam-public.ebi.ac.uk",
    user="rfamro",
    port ="4497",
    database="Rfam"
)
cursor = rfam_db.cursor()


print("Fetching ncRNA IDs and labels from Rfam")
#writing out what types of ncRNA we want here just so it will take up less space in the query string below
rna_types = ["rRNA", "snRNA", "sRNA", "tRNA", "ribozyme", "snoRNA", "miRNA", "lncRNA"]
type_condition = "(" + str.join(" or ", ["(f.type like \'%{x}%\')".format(x=x) for x in rna_types]) + ")"

query = """select distinct concat(fr.rfamseq_acc,'/',fr.seq_start,'-',fr.seq_end), f.type 
from full_region fr, family f 
where f.rfam_acc=fr.rfam_acc 
and fr.is_significant=1 
and {t};""".format(t=type_condition)

cursor.execute(query)

rows = cursor.fetchall()
fp = open('./ncrna.csv', 'w')
accessionFile = csv.writer(fp)
accessionFile.writerows(rows)
fp.close()

rfam_db.close()

print("Adding Sequences to dataset")
# Getting RNA sequences from Rfam.fa file
gene_info = pd.read_csv('./ncrna.csv', names=['accession', 'type'])

ids = set(gene_info['accession'].values)
seq_df_size = gene_info.shape[0]
seq_df = pd.DataFrame({"accession" : ["" for i in range(seq_df_size)], "sequence": ["" for i in range(seq_df_size)]})

with open('./Rfam.fa', 'r') as f:
    j = 0
    for line in f:
        #line beginning with '>' will be identifying information about sequence
        # First "word" in this line will be the same accession that
        # is in the gene info datafram
        if line[0] == '>':
            s = line.split()[0][1:]
        else:
            #Next line is sequence for the above accession
            seq = line
            # Only add sequence to dataframe if it's in the gene info dataframe
            if s in ids:
                seq_df.iloc[j, :] = [s,seq]
                ids.remove(s)
                j +=1
            if len(s) == 0:
                break

# Join sequence dataframe and gene info dataframe on accession
merged = pd.merge(seq_df, gene_info, left_on="accession", right_on="accession", how="inner")
merged = merged.rename(columns={"sequence_x": "sequence"})

merged.to_csv('./ncrna.csv', index=False)
