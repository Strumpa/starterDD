import os
import sys
import yaml
import math

current_dir = os.path.dirname(os.path.abspath(__file__))
root_dir = os.path.abspath(os.path.join(current_dir, '..'))
if root_dir not in sys.path:
    sys.path.insert(0, root_dir)

from starterDD.GeometryAnalysis.cartesian_geometry_analysis import CartesianGeometricAnalyser
from starterDD.InterfaceToDD.donjon_module_calls import DonjonTHM1DProcedure
# NOUVEL IMPORT : On importe le CoreModel
from starterDD.DDModel.DonjonModel import CoreModel 

def generer_nom_cas(pdrop, power_kw, type_profil):
    nom = "dfm"
    if pdrop == 1: nom += "p"
    nom += str(int(power_kw))
    if type_profil.lower() in ['cosinus', 'cos', 'c']: nom += "c"
    elif type_profil.lower() in ['sinus', 'sin', 's']: nom += "s"
    else: nom += "u"
    nom += "_v"
    return nom

def generer_profil_puissance(type_profil, nz):
    profil = []
    for i in range(nz):
        z_norm = (i + 0.5) / nz 
        if type_profil == 'cosinus':
            # Quart de cosinus : max en bas (z=0), nul en haut (z=1)
            val = math.cos((math.pi / 2.0) * z_norm) 
        elif type_profil == 'sinus':
            # Demi-sinus : cloche symétrique max au centre
            val = math.sin(math.pi * z_norm)
        else:
            val = 1.0
        profil.append(val)
    
    moyenne = sum(profil) / nz
    return [v / moyenne for v in profil]

def main():
    yaml_path = os.path.join(root_dir, 'data', 'BWRProgressionProblems', 'GE14_inputs', 'CORE', 'GEOM_4x4_mini_CORE.yaml')
    output_dir = os.path.join(current_dir, 'outputs', 'procs')
    
    if not os.path.exists(yaml_path):
        print(f"Erreur : Fichier {yaml_path} introuvable.")
        return

    # Reference values (PINLET, EPSOUT) for non-regression assertions.
    REFERENCE_VALS = {
        'dfmp10s_v': (7.23205200E+06, 6.43932670E-02),
        'dfmp20s_v': (7.23644300E+06, 4.47797000E-01),
        'dfmp40s_v': (7.25449950E+06, 6.77849114E-01),
        'dfmp10c_v': (7.23193850E+06, 6.07140325E-02),
        'dfmp20c_v': (7.23913300E+06, 4.56152081E-01),
        'dfmp40c_v': (7.26596800E+06, 6.83239937E-01),
        'dfm10c_v':  (7.20000000E+06, 6.04344234E-02),
        'dfm10s_v':  (7.20000000E+06, 6.40267283E-02),
        'dfm20c_v':  (7.20000000E+06, 4.43293720E-01),
        'dfm20s_v':  (7.20000000E+06, 4.34969038E-01),
        'dfm40c_v':  (7.20000000E+06, 6.60635710E-01),
        'dfm40s_v':  (7.20000000E+06, 6.55953765E-01),
    }

    print("--- Initialisation du modèle et de l'analyseur géométrique ---")
    core_pos = (1, 1)
    
    # --- MODIFICATIONS ICI ---
    # 1. Séparer le dossier et le nom du fichier
    path_to_configs = os.path.dirname(yaml_path)
    core_desc_file = os.path.basename(yaml_path)
    
    # 2. Instancier le modèle du cœur
    core_model = CoreModel(
        name="GE14_mini_core", 
        path_to_yaml_configs=path_to_configs, 
        core_description_yaml=core_desc_file
    )
    
    # 3. Créer les modèles d'assemblage en mémoire (Très important !)
    core_model.createAssemblyModels()
    
    # 4. On passe l'objet `core_model` à l'analyseur au lieu du chemin YAML
    analyser = CartesianGeometricAnalyser(core_model=core_model, core_i=core_pos[0], core_j=core_pos[1])
    # -------------------------

    # --- LA MATRICE DE TESTS ---
    puissances_a_tester = [0, 10.0, 20.0, 40.0]
    profils_a_tester = ['cosinus', 'sinus']
    pdrop_options = [0, 1]
    nz = 40
    dfm = 1

    print("--- Lancement de la génération par lots ---")
    
    # On boucle sur toutes les combinaisons possibles !
    for pdrop in pdrop_options:
        for power in puissances_a_tester:
            for profil_type in profils_a_tester:
                
                # 1. Génération du nom automatisé (ex: dfmp20c)
                case_name = generer_nom_cas(pdrop, power, profil_type)
                
                # 2. Génération du profil de puissance mathématique
                axial_pform = generer_profil_puissance(profil_type, nz)
                
                # 3. Création de la procédure
                procedure = DonjonTHM1DProcedure(
                    analyser=analyser,
                    nz=nz,
                    power_kw=power,
                    axial_pform=axial_pform, # On passe le vrai profil ici !
                    pdrop=pdrop,
                    dfm=dfm
                )
                
                # 4. Écriture du fichier
                refs = REFERENCE_VALS.get(case_name, (None, None))
                c2m_path = procedure.write_to_c2m(output_dir, case_name,
                                                   pinlet_ref=refs[0], epsout_ref=refs[1])
                print(f"✔️  Généré : {case_name}.c2m")

if __name__ == "__main__":
    main()