# Collection of functions used to help with geometry building, that do not require glow imports.
# R.Guasch
# Date : 09/02/2026 [created]

def computeSantamarinaradii(fuel_radius, gap_radius, clad_radius, gadolinium=False):
    """
    Helper to define fuel region radii for fuel pins
    A. Santamarina recommendations:
    - UOX pins: 50%, 80%, 95% and 100% volume fractions
    - Gd2O3 pins: 20%, 40%, 60%, 80%, 95% and 100% volume fractions
    """
    pin_radii = []
    if gap_radius is not None and gap_radius < fuel_radius:
        # in case expansion gap is inner most region, add it as the first radius and then compute the fuel radii based on the remaining fuel volume after subtracting the gap volume
        pin_radii.append(gap_radius)
    if not gadolinium:
        pin_radii.append([
            (0.5**0.5) * fuel_radius,
            (0.8**0.5) * fuel_radius,
            (0.95**0.5) * fuel_radius,
            fuel_radius
        ])
    else:
        pin_radii.append([
            (0.2**0.5) * fuel_radius,
            (0.4**0.5) * fuel_radius,
            (0.6**0.5) * fuel_radius,
            (0.8**0.5) * fuel_radius,
            (0.95**0.5) * fuel_radius,
            fuel_radius
        ])
    if gap_radius is not None and gap_radius > fuel_radius:
        # if gap radius is inner most region, add the clad radius as the last radius after the fuel radii
        pin_radii.append(gap_radius)    
    if clad_radius is not None:
        pin_radii.append(clad_radius)
    if gap_radius is not None and gap_radius > fuel_radius and clad_radius is None:
        ## print warning as this would be (fuel zones + gap zone) without clad. 
        print("Warning : gap radius is larger than fuel radius but no clad radius provided, this would lead to a pin with fuel zones and outer gap zone but no clad. Please check the input radii.")

    return pin_radii