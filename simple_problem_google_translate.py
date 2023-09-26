#! / usr / bin / env python
# - * - coding: utf-8 - * -
# sky mechanics I: a simple problem
# author: Franciszek Humieja
# version 1.1 (Jan. 26, 2017)

from sys import argv
import process_data as comp
import math as m
from scipy.optimize import fsolve
import find_kat

def load_orbital_elements (FileName):
dElementyOrbalne_linie = comp.read_data (FileName) #Read lines of text from the input file
dElementyOrbitalne_bez_kom comment = comp. delete_comments (dElementyOrbitalne_linie) #Delete comments
dElementyOrbitalne = comp.uporzadkuj_z_rozdz_tekst (dElementyOrbitalne_no_komkom) #Separates the data into two lists: the first with the names of the celestial bodies, and the second with their el. orbital (still as thongs)
return dElementsOrbital

def load_objectname (FileName):
objectNames = load_orbital_elements (FileName) [0]
return objectNames

def load_orbital_element_elements (FileName):
dElementyOrbalne = load_orbital_elements (FileName) [1]
n = flax (orbital elements)
elementsOrbital = []
for i in range (n):
elementsOrbital [i: i] = [[float (k) for k in dElementyOrbalne [i]]]
return itemsOrbital

def assign_watch_time (time):
"" "
Store the requested UTC time of observation as a list [y, m, d, g, m, s] where the second number is float and the rest are integer.
"" "
t_UTC_int = time [0: 5]
t_UTC_float = time [5]
t_UTC = [int (k) for k in t_UTC_int]
t_UTC [6: 6] = [t_UTC_float]
return t_UTC

def convert_deg_na_rad (kat):
kat_rad = kat / 180.0 * m.pi
return kat_rad

def convert_rad_na_deg (kat):
kat_deg = kat / m.pi * 180
return kat_deg

def convert_rad_na_h (kat):
cat_h = cat / m.pi * 12
return kat_h

def convert_time_to_JD (time):
"" "
Converts date and time from the Gregorian calendar to the Julian date. The argument of the function must be a list with the structure [year, month, day, hour, minute, second]. Source of the algorithm: English Wikipedia under the entry "Julian day".
"" "
year = time [0]
month = time [1]
day = time [2]
hour = time [3]
minute = time [4]
second = time [5]
a = int ((14 - month) / 12.) # January and February are the months of the previous year in JD; the variable is a 0-1 switch when these two months occur
y = year + 4800 - a # Year number from the beginning of the count
m = month + 12 * a - 3 # Month number from 0 (March) to 11 (February)
JDN = day + int ((153 * m + 2) / 5.) + 365 * y + int (y / 4.) - int (y / 100.) + Int (y / 400.) - 32045 # Second from the left component of the sum is taken from the different number of days in individual months
JD = JDN + (12 o'clock) / 24. + minute / 1440. + second / 86400.
return JD

def convert_JD_na_JED (JD):
"" "
Converts the JD (Julian Date) time, which is the universal UTC time, to the JED (Julian Ephemeris Date) time, which is the Earth time TT.
"" "
Dt = 69.184 # Difference in seconds between TT and UTC for 2017
JED = JD + Dt / 86400.
return JED

def load_lground_position (filename):
dPozosLand_linie = comp.read_data (filename)
dPozyZiemi_bez_kom comment = comp. deletion_comments (dPozyZiemi_linie)
dPozyZiemi = comp. order (dPozyZiemi_no_comments)
return dGround Positions

def ground_spolish_ interpolation (filename, t, coordinate_number):
"" "
Interpolates one selected barycentric coordinate of the position of the Earth.
"" "
dLandPositions = load_ground_positions (filename)
t_min = int (t) + 0.5 # Here I assume that the data in the Earth position file is set at 00:00 on each subsequent day!
t_max = int (t) + 1.5 # as above
coefficient = comp.find_wartosci_po_wierszach (dPosition of the Earth, [t_min, t_max], coordinate_number)
a = (coefficient [1] - coefficient [0]) / (t_max - t_min)
b = (coefficient [0] * t_max - coefficient [1] * t_min) / (t_max - t_min)
coordinate = a * t + b
return coordinate

def interpolate_ground_positions (filename, t):
"" "
Causes interpolation for all three barycentric Earth coordinates.
"" "
X_z = ground_spolordna_ interpoluj (filename, t, 1)
Y_z = ground_spolord_ interpoluj (filename, t, 2)
Z_z = interpoluj_wspolrzedna_ziemi (filename, t, 3)
return [X_z, Y_z, Z_z]

def compute_period (a):
d = 2 * m.pi / 0.01720209895 # Number of solar days in a sidereal year (a tropical year, due to the precession of the Earth, is shorter than a sidereal year by about 20.4 minutes)
T = m.sqrt (a ** 3) * d # Period in days; the formula is T = 2 * pi * sqrt (a ** 3 / (G * M_s)), but it is simplified when a is expressed in AU because G = 4 * pi ** 2 AU ** 3 yr ** (-2) M_s ** (- 1).
return T

def equation_on_anomalie_ eccentric (E, e, t, t_p, T, h):
equation = E - e * m.sin (E) - 2 * m.pi * (t - h - t_p) / T
return equation

def designate eccentric_anomalies (e, t, t_p, T, h):
"" "
Solves numerically the equation for the eccentric anomaly E. E should be in the range -pi to + pi.
"" "
E = fsolve (equation_na_anomalie_mimosrodowa, 0.0, (e, t, t_p, T, h))
return E [0]

def designate_anomalies_real (E, e):
phi = 2 * m.atan (m.sqrt ((1.0 + e) ​​/ (1.0 - e)) * m.tan (E / 2.0)) # Simple arctg formula since phi finds ranges from -pi to + pi
return phi

def designate_odleglosc_od_barycentrum (a, e, E):
R = a * (1 - e * m.cos (E))
return R

def designate_polozenie_barycentyczne (R, Omega, phi, omega, i):
X = R * (m.cos (Omega) * m.cos (phi + omega) - m.sin (Omega) * m.cos (i) * m.sin (phi + omega))
Y = R * (m.sin (Omega) * m.cos (phi + omega) + m.cos (Omega) * m.cos (i) * m.sin (phi + omega))
Z = R * m.sin (i) * m.sin (phi + omega)
return [X, Y, Z]

def designate_geocentric_position (X, X_z):
epsilon_deg = 23.439292 # slope of the ecliptic to the equator [deg]
epsilon = convert_deg_na_rad (epsilon_deg) # [rad]
x = X [0] - X_z [0]
y = (X [1] - X_z [1]) * m.cos (epsilon) - (X [2] - X_z [2]) * m.sin (epsilon)
z = (X [1] - X_z [1]) * m.sin (epsilon) + (X [2] - X_z [2]) * m.cos (epsilon)
return [x, y, z]

def designate_distance_from_ground (x):
r = m.sqrt (x [0] ** 2 + x [1] ** 2 + x [2] ** 2)
return

def designate_delay_signal (r):
c_ms = 299792458.0 # speed of light in vacuum [m / s]
d = 86400.0 # converter s to d [s / d]
a = 149597870700.0 # m to AU converter [m / AU]
c = c_ms * d / a # [AU / d]
h = r / c
return h

def designate_rektascences_and_declinations (x, r):
delta_rad = m.asin (x [2] / r) #The declination varies from -pi / 2 to + pi / 2, so you can use the normal inverse sine function; [glad]
cos_alpha = x [0] / (r * m.cos (delta_rad))
sin_alpha = x [1] / (r * m.cos (delta_rad))
alpha_rad = find_kat.z_funkcji_tryg (sin_alpha, cos_alpha) #Rectascension varies between 0 and 2 * pi, so you can't use the normal inverse function; [glad]
alpha = convert_rad_na_h (alpha_rad) # [arch]
delta = convert_rad_na_deg (delta_rad) # [deg]
return [alpha, delta]

def save_ephemerides (Filename, ObjectNames, ObservationTimes, DistanceFrom Earth, Equivalent Coordinates):
dZapis = comp.zlacz_dane_z_jednej_linii (names of objects, observation times, distances from the ground, co-ordinates)
dSave [0: 0] = [["#name", "UTC (r", "m", "d", "g", "m", "s)", "distance [AU]", "RA [arch] "," DEC [deg] "]]
save_data (file name, dSave)

if __name__ == "__main__":
docl = 12 # Accuracy with which the values ​​of successive determined eccentric anomalies E are to be compared
fElementyOrbalne = argv [1] #Fetching the name of the input file with elements, which is an argument of calling the program
fEarth position = "earth_position.dat" # Name of the input file with the Earth's coordinates in the barycentric system
fEfemerydy = "ephemerides.dat" #Output file name with ephemeris
ObjectNames = load_Object_name (fElementsOrbital)
orbital_elements = load_wartosci_elementow_orbitalnych (fElementyOrbitalnych)
timesObservation = []
odleglosciOdZiemi = []
co-ordinates = []
n = len (orbital elements)
for k in range (n): # Loop through list items
e = Orbital [k] [0] # eccentricity
a = Orbital elements [k] [1] # major driveshaft [AU]
t_p = orbital elements [k] [2] # time of passage through the periapsis [JD]
Omega = recalculate_deg_na_rad (elementsOrbital [k] [3]) # length of the ascending node [rad]
i = recalculate_deg_na_rad (elementsOrbital [k] [4]) #inclination [rad]
omega = recalculate_deg_na_rad (elementsOrbital [k] [5]) #argument pericenter [rad]
t_UTC = assign_observation_time (elementsOrbital [k] [6:12]) # desired UTC time of observation [r, m, d, g, m, s]
t_JD = convert_time_to_JD (t_UTC) # [d]
t = convert_JD_to_JED (t_JD) # [d]
X_z = interpolate_ground_position (fEarth_position, t) # [AU, AU, AU]
T = calculate_period (a) # [d]
h = 0.0 # time of the light signal from the celestial body to the Earth --- correction to the formula for an eccentric anomaly
E_kontrolna = 0.0 # eccentric anomaly
end = False
while not end:
E = determine eccentric_anomalies (e, t, t_p, T, h) # [rad]
phi = designate_anomalie_real (E, e) # [rad]
R = determine_odleglosc_od_barycentrum (a, e, E) # [AU]
X = designate_polozenie_barycentyczne (R, Omega, phi, omega, i) # [AU, AU, AU]
x = designate_geocentric_position (X, X_z) # [AU, AU, AU]
r = designate_distance_from_ground (x) # [AU]
h = signal_delay (r) # [d]
if round (E - E_kontrolna, dokl) == 0.0: #Checking if the corrections affect the value of E
end = True
E_kontrolna = E
co-ordinatesRownikowe_k = designate_rektascensje_i_declinacje (x, r) # [arch, deg]
timesObserwacji [k: k] = [t_UTC]
odleglosciOdZiemi [k: k] = [r]
co-ordinate [k: k] = [co-ordinate_k]
save_ephemerides (fEfemerides, names of objects, times of observation, distances from the ground, co-ordinates)
