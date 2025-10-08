import matplotlib.pyplot as plt
import numpy as np

distrib_angles = [ 0.052, 0.157, 0.259, 0.358, 0.455, 0.544, 0.629, 0.709, 0.777, 0.841, 0.89, 0.934, 0.968, 0.993, 1, 0.997, 0.989, 0.967, 0.932, 0.892, 0.842, 0.776, 0.709, 0.628, 0.543, 0.456, 0.358, 0.256, 0.156,0.0529]


plt.figure(figsize=(10, 6))

num_bins = len(distrib_angles)
x = [(a+0.5)/num_bins*np.pi/2.0 for a in range(num_bins)]

print (x)
plt.plot(x, distrib_angles, 'x',label='Angular distribution')

# Assuming uniform bin width and the distribution is over 0 to pi/2
bin_centers = np.linspace(0, np.pi/2, num_bins)
cosine_fit = np.sin(2*bin_centers)

# Normalize the cosine fit to match the scale of the bar plot if necessary
# This simple normalization might not be accurate, a proper fit would involve optimization
cosine_fit_normalized = cosine_fit

plt.plot(bin_centers, cosine_fit_normalized, color='red', linestyle='--', label='Sine Fit')


plt.xlabel("angle (rad)")
plt.ylabel("Probability Density")
plt.title("Bar Plot of distrib_angles with Cosine Fit")
plt.legend()
plt.grid(axis='y')
plt.savefig("cosine_distrib.png",bbox_inches='tight')
