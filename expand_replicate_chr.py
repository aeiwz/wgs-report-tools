import pandas as pd

class ExpandReplicateChr:
    def __init__(self, df):
        self.df = df
    
    def expand_rows(self):
        # Ensure 'CHR_POS' is treated as a string and expand it
        self.df['CHR_POS'] = self.df['CHR_POS'].astype(str).str.replace(' x ', ';')
        expanded_df = self.df.assign(CHR_POS=self.df['CHR_POS'].str.split(';')).explode('CHR_POS')

        return expanded_df

    def process_chromosome_ids(self):
        # Ensure 'CHR_ID' is a string and handle ';' separator
        self.df['CHR_ID'] = self.df['CHR_ID'].astype(str).str.split(';').str[0]

        # Ensure consistent formatting for 'CHR_ID'
        self.df['CHR_ID'] = self.df['CHR_ID'].str.split(' x ').str[0]
        self.df['CHR_ID'] = 'chr' + self.df['CHR_ID'].str.replace('.0', '', regex=False)

        return self.df

    def preprocess_gwas_data(self):
        print ('Missing values have been found:', self.df['CHR_POS'].isnull().sum())
        self.df.dropna(subset=['CHR_POS'], inplace=True)
        print ('Missing values were droped')
        self.df = self.expand_rows()
        print ('Rows were expanded')
        self.df = self.process_chromosome_ids()
        print ('Chromosome IDs were processed')

        # Convert 'CHR_POS' to numeric, coercing errors to NaN, then filling with 0 and converting to int
        self.df['CHR_POS'] = pd.to_numeric(self.df['CHR_POS'], errors='coerce').fillna(0).astype(int)

        return self.df
    
if __name__ == '__main__':
    # Load and process the dataset
    df = pd.read_csv('data/gwas_database_with_description.csv.gz', compression='gzip', low_memory=False)
    expand_replicate_chr = ExpandReplicateChr(df)
    df2 = expand_replicate_chr.preprocess_gwas_data()
    df2.to_csv('data/gwas_database_with_description_expanded.csv.gz', compression='gzip', index=False)

