from magnetic_field_calculator import MagneticFieldCalculator
import numpy as np

geographic_times=np.loadtxt('latitude_vs_longitude_times.tsv', delimiter="\t", dtype=str, usecols=range(361))
geographic_longitudes=geographic_times[0,1:]
geographic_latitudes=geographic_times[1:,0]
calculator = MagneticFieldCalculator()

file=np.loadtxt("latitude_distribution_oneway.txt", delimiter="\t", dtype=str)
orientations=file[:,2]
print(orientations)
results=[]

distribution_file3=open("geomagnetic_headings.txt", "w")

for i in range(len(geographic_latitudes)):
    for j in range(len(geographic_longitudes)):
        fraction=geographic_times[i+1,j+1]
        # print(j)
        # print(geographic_latitudes[i], geographic_longitudes[j], fraction)
        # print(geographic_latitudes[i], geographic_longitudes[j], fraction)
        calc=calculator.calculate(latitude=geographic_latitudes[i], longitude=geographic_longitudes[j], altitude=450, date="2026-07-14")
        field_value = calc['field-value']
        declination = field_value['declination']

        print(geographic_latitudes[i], float(geographic_longitudes[j]), float(fraction), orientations[i], declination['value'])
        distribution_file3.write(f"{geographic_latitudes[i]}\t{float(geographic_longitudes[j])}\t{float(fraction)}\t{orientations[i]}\t{declination['value']}\n")

distribution_file3.close()