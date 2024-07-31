#!/usr/bin/env python3

import subprocess
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

# Define the command
command = "./abcore 1 ../bidata/rmwiki 4 0 0 | grep estis"

# Initialize lists to store results
naive_estimates = []
bs_estimates = []
advss = []
advds = []

# Run the command 100 times
for _ in range(1000):
    # Execute the command and capture the output
    output = subprocess.check_output(command, shell=True).decode("utf-8")

    # Split the output lines
    lines = output.strip().split("\n")

    # Extract the estimates from the output
    cnter = 0
    for line in lines:
        if "naive estis" in line:
            naive_estimate = int(line.split(":")[-1].strip())
            naive_estimates.append(naive_estimate)
        elif "bs estis" in line:
            bs_estimate = float(line.split(":")[-1].strip())
            bs_estimates.append(bs_estimate)
        elif "adv estis" in line:
            adv_estimate = float(line.split(":")[-1].strip())
            print("adv estis: ", adv_estimate)
            if cnter ==0:
                advss.append(adv_estimate)
            if cnter ==2:
                advds.append(adv_estimate)
            cnter = cnter +1

# Convert lists to numpy arrays
naive_estimates = np.array(naive_estimates)
bs_estimates = np.array(bs_estimates)
advss = np.array(advss)
advds = np.array(advds)

print(naive_estimates)
print(bs_estimates)
print(advss)
print(advds)

# Create probability density curves
sns.set_style("whitegrid")
plt.figure(figsize=(10, 6))
sns.kdeplot(naive_estimates, color='blue', label='Naive', shade=True)
sns.kdeplot(bs_estimates, color='orange', label='OneR', shade=True)
sns.kdeplot(advss, color='green', label='MultiR-SS', shade=True)
sns.kdeplot(advds, color='red', label='MultiR-DS', shade=True)


# Add vertical lines for mean values
# Calculate mean and standard deviation for each array
# naive_mean = np.mean(naive_estimates)
# bs_mean = np.mean(bs_estimates)
# first_adv_mean = np.mean(advss)
# last_adv_mean = np.mean(advds)
# plt.axvline(x=naive_mean, color='blue', linestyle='--', label='Naive Mean')
# plt.axvline(x=bs_mean, color='orange', linestyle='--', label='OneR Mean')
# plt.axvline(x=first_adv_mean, color='green', linestyle='--', label='MultiR-SS Mean')
# plt.axvline(x=last_adv_mean, color='red', linestyle='--', label='MultiR-DS Mean')

# Add a vertical line for the actual count of 9
plt.axvline(x=2, color='black', linestyle='--', linewidth=1, label='True Count')



plt.xlabel('Estimate')
plt.ylabel('Probability Density')
plt.title('Probability Density of Estimates')
plt.legend()
plt.grid(True)
# plt.show()


plt.savefig('distribution.pdf', format='pdf')

plt.close()
# plt.show()
