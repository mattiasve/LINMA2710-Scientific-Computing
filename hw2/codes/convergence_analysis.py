#%%
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt 
import os 
wd = os.getcwd()

df0 = pd.read_csv(wd+'/convergence_analysis_sol/sol_conan0.txt', sep=" ", header=None)
df1 = pd.read_csv(wd+'/convergence_analysis_sol/sol_conan1.txt', sep=" ", header=None) 
df2 = pd.read_csv(wd+'/convergence_analysis_sol/sol_conan2.txt', sep=" ", header=None) 
df3 = pd.read_csv(wd+'/convergence_analysis_sol/sol_conan3.txt', sep=" ", header=None) 
df4 = pd.read_csv(wd+'/convergence_analysis_sol/sol_conan4.txt', sep=" ", header=None) 
# df5 = pd.read_csv(wd+'/convergence_analysis_sol/sol_conan5.txt', sep=" ", header=None) 
# df6 = pd.read_csv(wd+'/convergence_analysis_sol/sol_conan6.txt', sep=" ", header=None) 
# df7 = pd.read_csv(wd+'/convergence_analysis_sol/sol_conan7.txt', sep=" ", header=None) 
# df8 = pd.read_csv(wd+'/convergence_analysis_sol/sol_conan8.txt', sep=" ", header=None) 
# df9 = pd.read_csv(wd+'/convergence_analysis_sol/sol_conan9.txt', sep=" ", header=None) 

fig = plt.figure(figsize=(8,8))
plt.scatter(df4[0], df4[1], s=10, label='res = 0.25')
plt.scatter(df3[0], df3[1], s=10, label='res = 0.5')
plt.scatter(df2[0], df2[1], s=10, label='res = 1')
plt.scatter(df1[0], df1[1], s=10, label='res = 2')
plt.scatter(df0[0], df0[1], s=10, label='ref. points | res = 4')
plt.legend(ncol=5, loc='upper center', bbox_to_anchor=(0.5, -0.05))

plt.savefig('./plots/mesh_resolution.svg', bbox_inches='tight')
plt.savefig('./plots/mesh_resolution.pdf', bbox_inches='tight', dpi=300)
plt.show()
#%%
# metrics computation
ref_pts = df0.iloc[:, :2]
ref_df1 = pd.merge(ref_pts, df1, on=[0, 1])
ref_df2 = pd.merge(ref_pts, df2, on=[0, 1])
ref_df3 = pd.merge(ref_pts, df3, on=[0, 1])
ref_df4 = pd.merge(ref_pts, df4, on=[0, 1])

errors = []
diff_nb_pts = []
lst_ref = [df0, ref_df1, ref_df2, ref_df3, ref_df4]
nb_points = [1022, 1070, 1306, 2246, 6006]

for i in range(len(lst_ref)-1):
    dfi = lst_ref[i]
    dfiplus1 = lst_ref[i+1]
    sum_error = 0
    for j in range(len(ref_pts)):
        sum_error += (dfiplus1.iloc[j][2] - dfi.iloc[j][2])**2
    rmse = np.sqrt(sum_error/len(ref_pts))
    errors.append(rmse)

for i in range(len(nb_points)-1):
    diff_nb_pts.append(float(nb_points[i+1]-nb_points[i]))

print(diff_nb_pts)
dummy_pts = np.reciprocal(diff_nb_pts)

# plots
plt.figure(dpi=300)
plt.plot(diff_nb_pts, errors, '-o', label='data')
#plt.plot(diff_nb_pts, dummy_pts, '-o', label=r'$y=1/x$')
plt.xlabel('Additional number of points for two successive mesh resolution')
plt.ylabel('RMSE [-]')
plt.grid(linestyle=':')
plt.savefig('./plots/convergence_analysis.svg', bbox_inches='tight')
plt.savefig('./plots/convergence_analysis.pdf', bbox_inches='tight', dpi=300)
plt.show()

plt.figure(dpi=300)
plt.loglog(diff_nb_pts, errors, '-o', label='data')
plt.loglog(diff_nb_pts, dummy_pts, '-o', label=r'$y=1/x$')
plt.xlabel('Log additional number of points for two successive mesh resolution')
plt.ylabel('Log RMSE [-]')
plt.grid(linestyle=':')
plt.legend(loc='best')
plt.savefig('./plots/convergence_analysis_log.svg', bbox_inches='tight')
plt.savefig('./plots/convergence_analysis_log.pdf', bbox_inches='tight', dpi=300)
plt.show()

# %%
