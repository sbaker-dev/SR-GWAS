from miscSupports import terminal_time
import statsmodels.api as sm
from pathlib import Path
import pandas as pd

# Setup working directory
working_directory = Path(Path().resolve(), "example_notebook.ipynb")
if not working_directory.exists():
    raise ReferenceError("WARNING: PATH TO WORKING DIRECTORY COULD NOT BE ESTABLISHED")
else:
    working_directory = working_directory.parent

# Covariant list
covariant_list = ["Gender", "Age", "PC1", "PC2", "PC3", "PC4", "PC5", "PC6", "PC7", "PC8", "PC9", "PC10", "Constant"]

# Setup the database
database = pd.read_csv(Path(working_directory, "Data", "CovariantSnp.csv"))
database["Constant"] = [1 for i in range(len(database))]
print(database)
print(f"Loaded Environment {terminal_time()}")

# Model 1
# regress BMI on G, sex, YoB, PCs
model = sm.OLS(database["BMI"], database[["rs012"] + covariant_list], missing="drop").fit()
print(model.summary())

# Model 2
# Residualise BMI and then regress residualised BMI on G

res = sm.OLS(database["BMI"], database[covariant_list], missing="drop").fit()
print(res.resid)

model = sm.OLS(res.resid, database[["rs012", "Constant"]]).fit()
print(model.summary())


# Model 3
# Residualise G and then regress BMI on residualised G
g_res = sm.OLS(database["rs012"], database[covariant_list + ["Constant"]], missing="drop").fit()
g_res = pd.concat([pd.DataFrame(g_res.resid, columns=["rs012"]), database["Constant"]], axis=1)
print(g_res)

model = sm.OLS(database["BMI"], g_res, missing="drop").fit()
print(model.summary())

# Model 4
# Regress of residualised BMI on residualised G
model = sm.OLS(res.resid, g_res, missing="drop").fit()
print(model.summary())
