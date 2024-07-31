
import uproot
import numpy as np
import math
import ROOT

input_file_path = "/eos/user/j/jaeyserm/public/mumuH_rec_16478_29.root"
file = uproot.open(input_file_path)

events = file["events"]
vtx_barrel = events[209]

radii = [14, 36, 58]

x_data = vtx_barrel["VertexBarrelCollection.position.x"].array()
y_data = vtx_barrel["VertexBarrelCollection.position.y"].array()
z_data = vtx_barrel["VertexBarrelCollection.position.z"].array()

px_data = vtx_barrel["VertexBarrelCollection.momentum.x"].array()
py_data = vtx_barrel["VertexBarrelCollection.momentum.y"].array()
pz_data = vtx_barrel["VertexBarrelCollection.momentum.z"].array()

cellid_data = vtx_barrel["VertexBarrelCollection.cellID"].array()
edep_data = vtx_barrel["VertexBarrelCollection.EDep"].array()
time_data = vtx_barrel["VertexBarrelCollection.time"].array()
pathlength_data = vtx_barrel["VertexBarrelCollection.pathLength"].array()

def radius(r):
    for true_r in radii:
        if abs(r - true_r) < 5:
            return true_r
    raise ValueError(f"Not close enough to any of the layers {np.sqrt(true_r)}")

def phi(x,y):
    """
    Calculates phi of particle.
    Inputs: x,y floats.
    Output: phi, float representing angle in radians from 0 to 2 pi.
    """
    phi = math.atan(y/x)
    if x < 0:
        phi +=  math.pi
    elif y < 0:
        phi += 2*math.pi
    return phi

def theta(x,y,z):
    """
    Calculates theta of particle.
    Inputs: x,y,z floats.
    Output: theta, float representing angle in radians from 0 to pi.
    """
    return math.acos(z/np.sqrt(x**2 + y**2 + z**2))

hit_map = {r: {cos: {azimuthal: [0]*1000 for azimuthal in range(0, 360, 3)} for cos in range(-100, 100, 2)} for r in radii}
for i in range(1000):

    hits = {coord: [] for coord in radii}
    event_x = x_data[i]
    event_y = y_data[i]
    event_z = z_data[i]

    for j in range(len(event_x)):
        x = event_x[j]
        y = event_y[j]
        r = np.sqrt(x**2 + y**2)
        hits[radius(r)].append(j)

    for r in radii:
        for index in hits[r]:
            x = event_x[index]
            y = event_y[index]
            z = event_z[index]
            azimuthal = phi(x, y) * (180 / math.pi)
            polar = theta(x, y, z)
            cos = 100 * math.cos(polar)
            hit_map[r][int((cos // 2) * 2)][int((azimuthal // 3) * 3)][i] += 1

for layer_index in range(3):

    r = radii[layer_index]
    hist = ROOT.TH2F("hit map", f"Guinea Pig Layer {layer_index + 1} Module Hits", \
                    100, -100, 100, 120, 0, 360)
    hist.SetTitle(f"Layer {layer_index + 1} Hits per Bunch Crossing;Cosine Theta;Azimuthal Angle (deg)")

    for cos in range(-100, 100, 2):
        for azimuthal in range(0, 360, 3):
            hits = np.mean(hit_map[r][cos][azimuthal])
            hist.SetBinContent(((cos + 100) // 2) + 1, (azimuthal // 3) + 1, hits)

    hist.SetStats(0)
    canvas = ROOT.TCanvas("hit map", f"Layer {layer_index + 1} Hits")
    hist.Draw("colz")
    canvas.Update()
    canvas.SaveAs(f"../plots/hit_rates/higgs/higgs_layer{layer_index + 1}_hit_rate_test.png")
