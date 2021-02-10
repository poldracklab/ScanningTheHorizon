
# %%
import seaborn as sns
import pandas as pd
import numpy as np

# %%
df = pd.read_csv('group_n_data.csv')

# %%
df['log_n'] = np.log(df['n'])
sns.scatterplot(x='year', y='log_n', data=df)
sns.lineplot(x='year', y='log_n', data=df,
             ci=None, estimator="median")
# %%
medians = df.groupby('year').median()
# %%
