import argparse
import os
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
from starterDD.GeometryAnalysis.geometry_analysis import GeometricAnalyser


RACINE_PROJET = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
DOSSIER_OUTPUTS = os.path.join(RACINE_PROJET, "tests", "outputs")
 
if __name__ == "__main__":
    # Création du dossier outputs s'il n'existe pas encore
    os.makedirs(DOSSIER_OUTPUTS, exist_ok=True)
    
    parser = argparse.ArgumentParser(description="Test de l'Analyseur Géométrique via fichier CORE")
    
    # 1. Le fichier CORE est obligatoire en premier
    parser.add_argument("core_file", help="Chemin vers le fichier CORE YAML à tester")
    
    # 2. Position de l'assemblage dans le cœur (Défaut: 1,1)
    parser.add_argument("--core_pos", nargs=2, type=int, default=[1, 1], metavar=('CORE_I', 'CORE_J'), help="Position de l'assemblage dans le coeur")
    
    parser.add_argument("--x1", type=float, default=None)
    parser.add_argument("--x2", type=float, default=None)
    parser.add_argument("--y1", type=float, default=None)
    parser.add_argument("--y2", type=float, default=None)
    parser.add_argument("--z1", type=float, default=None)
    parser.add_argument("--z2", type=float, default=None)
    
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("--cv_z", action="store_true")
    group.add_argument("--cv_x", action="store_true")
    group.add_argument("--cv_y", action="store_true")
    group.add_argument("--rod_z", nargs=2, type=int, metavar=('I', 'J'))
    group.add_argument("--rod_x", nargs=2, type=int, metavar=('I', 'J'))
    group.add_argument("--rod_y", nargs=2, type=int, metavar=('I', 'J'))
    group.add_argument("--water_z", nargs=2, type=int, metavar=('I', 'J'))
    group.add_argument("--water_x", nargs=2, type=int, metavar=('I', 'J'))
    group.add_argument("--water_y", nargs=2, type=int, metavar=('I', 'J'))
    group.add_argument("--mesh_z", choices=['fuel', 'water', 'regular'])
    group.add_argument("--mesh_x", choices=['fuel', 'water', 'regular'])
    group.add_argument("--mesh_y", choices=['fuel', 'water', 'regular'])
    
    parser.add_argument("--n", type=int, default=10)
    parser.add_argument("--n_z", type=int, default=10)
    parser.add_argument("--plot", action="store_true")
    
    parser.add_argument("--plot_z", action="store_true")
    parser.add_argument("--plot_x", action="store_true")
    parser.add_argument("--plot_y", action="store_true")
    parser.add_argument("--h", type=float, default=10.0)
    parser.add_argument("--p", type=float, default=1.0)

    
    args = parser.parse_args()
    
    # On initialise l'API avec le fichier CORE et la position demandée
    analyser = GeometricAnalyser(args.core_file, core_i=args.core_pos[0], core_j=args.core_pos[1])
    x_min, x_max = analyser.get_x_global_bounds()
    y_min, y_max = analyser.get_y_global_bounds()
    z_min, z_max = analyser.get_z_global_bounds()
    x1 = args.x1 if args.x1 is not None else x_min
    x2 = args.x2 if args.x2 is not None else x_max
    y1 = args.y1 if args.y1 is not None else y_min
    y2 = args.y2 if args.y2 is not None else y_max
    z1 = args.z1 if args.z1 is not None else z_min
    z2 = args.z2 if args.z2 is not None else z_max
    section_type = ['cv',[x1,y1,x2,y2]] if args.cv_z else ['cv', [y1,z1,y2,z2]] if args.cv_x else ['cv', [x1,z1,x2,z2]] if args.cv_y else ['rod', args.rod_z] if args.rod_z else ['rod', args.rod_x] if args.rod_x else ['rod', args.rod_y] if args.rod_y else ['water', args.water_z] if args.water_z else ['water', args.water_x] if args.water_x else ['water', args.water_y] if args.water_y else [None, None]
    mesh_type_z= ['regular', args.n] if args.mesh_z == 'regular' else [args.mesh_z, None]
    mesh_type_x= ['regular', args.n] if args.mesh_x == 'regular' else [args.mesh_x, None]
    mesh_type_y= ['regular', args.n] if args.mesh_y == 'regular' else [args.mesh_y, None]
    n_z = args.n_z
    h= args.h
    p= args.p

    if args.plot_z:
        if not (args.cv_z or args.rod_z or args.water_z):
            print("Erreur : Spécifiez une cible (--cv_z, --rod_z, ou --water_z) pour --plot_z.")
            exit(1)
        z_coords, porosities, a_cools, dhs, phs, cible_str= analyser.execute_profile_z(section_type, h, p, z_min, z_max)
        fig, axs = plt.subplots(4, 1, figsize=(10, 12), sharex=True)
        axs[0].plot(z_coords, porosities, label=r'Porosité ($\phi$)', color='blue')
        axs[0].set_ylabel(r"Porosité $\phi$ (-)")
        axs[1].plot(z_coords, a_cools, label=r'Acool (cm²)', color='orange')
        axs[1].set_ylabel(r"Surface de passage fluide $A_{cool}$ (cm²)")
        axs[2].plot(z_coords, dhs, label='Dh', color='green')
        axs[2].set_ylabel("Diamètre Hydraulique Dh (cm)")
        axs[3].plot(z_coords, phs, label='Ph', color='red')
        axs[3].set_ylabel("Périmètre de Chauffe Ph (cm)")
        axs[3].set_xlabel("Altitude axiale Z (cm)")
        
        titre = f"Profil Axial (Fenêtre h={args.h}cm, Pas p={args.p}cm) - {cible_str} [Assm: {args.core_pos[0]},{args.core_pos[1]}]"
        
        for ax in axs: 
            ax.grid(True)
            ax.legend()
            
        plt.suptitle(titre)
        plt.tight_layout()
        
        # SAUVEGARDE AU LIEU D'AFFICHER
        nom_fichier = f"profil_z_{cible_str}_assm{args.core_pos[0]}{args.core_pos[1]}_h{args.h}_p{args.p}.png"
        chemin_complet = os.path.join(DOSSIER_OUTPUTS, nom_fichier)
        plt.savefig(chemin_complet, dpi=300, bbox_inches='tight')
        plt.close() 
        print(f"Graphique sauvegardé avec succès dans :\n -> {chemin_complet}\n")
    
    elif args.plot_x:
        if not (args.cv_x or args.rod_x or args.water_x):
            print("Erreur : Spécifiez une cible (--cv_x, --rod_x, ou --water_x) pour --plot_x.")
            exit(1)
        x_coords, porosities, dhs, cible_str = analyser.execute_profile_x(section_type, z1, z2, p)
        fig, axs = plt.subplots(2, 1, figsize=(10, 12), sharex=True)
        axs[0].plot(x_coords, porosities, label=r'Porosité ($\phi$)', color='blue')
        axs[0].set_ylabel(r"Porosité $\phi$ (-)")
        axs[1].plot(x_coords, dhs, label='Dh', color='green')
        axs[1].set_ylabel("Diamètre Hydraulique Dh (cm)")
        axs[1].set_xlabel("Coordonnée X (cm)")
        
        titre = f"Profil Horizontal (Pas p={args.p}cm) - {cible_str} [Assm: {args.core_pos[0]},{args.core_pos[1]}]"
        
        for ax in axs: 
            ax.grid(True)
            ax.legend()
            
        plt.suptitle(titre)
        plt.tight_layout()
        
        # SAUVEGARDE AU LIEU D'AFFICHER
        nom_fichier = f"profil_x_{cible_str}_Z{z1}_{z2}_assm{args.core_pos[0]}{args.core_pos[1]}_p{args.p}.png"
        chemin_complet = os.path.join(DOSSIER_OUTPUTS, nom_fichier)
        plt.savefig(chemin_complet, dpi=300, bbox_inches='tight')
        plt.close() 
        print(f"Graphique sauvegardé avec succès dans :\n -> {chemin_complet}\n")
    
    elif args.plot_y:
        if not (args.cv_y or args.rod_y or args.water_y):
            print("Erreur : Spécifiez une cible (--cv_y, --rod_y, ou --water_y) pour --plot_y.")
            exit(1)
        y_coords, porosities, dhs, cible_str = analyser.execute_profile_y(section_type, z1, z2, p)
        fig, axs = plt.subplots(2, 1, figsize=(10, 12), sharex=True)
        axs[0].plot(y_coords, porosities, label=r'Porosité ($\phi$)', color='blue')
        axs[0].set_ylabel(r"Porosité $\phi$ (-)")
        axs[1].plot(y_coords, dhs, label='Dh', color='green')
        axs[1].set_ylabel("Diamètre Hydraulique Dh (cm)")
        axs[1].set_xlabel("Coordonnée Y (cm)")
        
        titre = f"Profil Horizontal (Pas p={args.p}cm) - {cible_str} [Assm: {args.core_pos[0]},{args.core_pos[1]}]"

        for ax in axs: 
            ax.grid(True)
            ax.legend()
        
        plt.suptitle(titre)
        plt.tight_layout()

        # SAUVEGARDE AU LIEU D'AFFICHER
        nom_fichier = f"profil_horizontal_{cible_str}_Z{z1}_{z2}_assm{args.core_pos[0]}{args.core_pos[1]}_p{args.p}.png"
        chemin_complet = os.path.join(DOSSIER_OUTPUTS, nom_fichier)
        plt.savefig(chemin_complet, dpi=300, bbox_inches='tight')
        plt.close() 
        print(f"Graphique sauvegardé avec succès dans :\n -> {chemin_complet}\n")

    elif args.mesh_z:
        mat_p, mat_dh, mat_ph = analyser.execute_mesh_z(mesh_type_z, z1, z2)
        print("--- Matrice de Porosité (phi) ---")
        print(np.round(np.array(mat_p), 4)) # Arrondi à 4 décimales pour la lisibilité
        
        print("\n--- Matrice de Diamètre Hydraulique (Dh) ---")
        print(np.round(np.array(mat_dh), 4))
        
        print("\n--- Matrice de Périmètre Chauffant (Ph) ---")
        print(np.round(np.array(mat_ph), 4))
        if args.plot:
            matrice_numpy = np.array(mat_p)
            plt.figure(figsize=(10, 8))
            matrice_a_dessiner = matrice_numpy[::-1]
            afficher_nombres = (len(mat_p[0]) <= 15)

            ass_geo = analyser.data_ref['ASSEMBLY_GEOMETRY']
            W_start = ass_geo['gap_wide'] + ass_geo['channel_box_thickness']
            W_end = ass_geo['assembly_pitch'] - W_start
            
            y_coords = np.linspace(W_end, W_start, len(matrice_numpy))
            step_y = max(1, len(y_coords) // 15) 
            labels_y = [f"{val:.2f}" if i % step_y == 0 else "" for i, val in enumerate(y_coords)]
            
            x_coords_phys = np.linspace(W_start, W_end, len(matrice_numpy[0]))
            step_x = max(1, len(x_coords_phys) // 10)
            labels_x = [f"{val:.2f}" if i % step_x == 0 else "" for i, val in enumerate(x_coords_phys)]
            
            sns.heatmap(matrice_a_dessiner, annot=afficher_nombres, fmt=".2f", 
                        cmap="viridis", cbar_kws={'label': r'Porosité $\phi$ (-)'},
                        yticklabels=labels_y, xticklabels=labels_x)
            
            plt.yticks(rotation=0) 
            plt.xticks(rotation=45) 
            
            plt.title(f"Carte de Porosité 3D (Entre Z={z1}cm et Z={z2}cm) - {args.mesh_z.capitalize()} - Assm: {args.core_pos[0]},{args.core_pos[1]}")
            plt.xlabel("Axe X (en cm)")
            plt.ylabel("Axe Y (en cm)")
            
            nom_fichier = f"heatmap_mesh_z_{args.mesh_z}_assm{args.core_pos[0]}{args.core_pos[1]}_Z{z1}_Z{z2}.png"
            chemin_complet = os.path.join(DOSSIER_OUTPUTS, nom_fichier)
            plt.savefig(chemin_complet, dpi=300, bbox_inches='tight')
            plt.close() 
            
            print(f"Heatmap sauvegardée avec succès dans :\n -> {chemin_complet}\n")
        
    elif args.mesh_x:
        mat_p, mat_dh = analyser.execute_mesh_x(mesh_type_x, x1, z1, z2, n_z)
        print("--- Matrice de Porosité (phi) ---")
        print(np.round(np.array(mat_p), 4)) # Arrondi à 4 décimales pour la lisibilité
        
        print("\n--- Matrice de Diamètre Hydraulique (Dh) ---")
        print(np.round(np.array(mat_dh), 4))

        if args.plot:
            matrice_numpy = np.array(mat_p)
            plt.figure(figsize=(10, 8))
            matrice_a_dessiner = matrice_numpy[::-1]
            afficher_nombres = (len(mat_p[0]) <= 15 and len(mat_p) <= 15)   

            ass_geo = analyser.data_ref['ASSEMBLY_GEOMETRY']
            W_start = ass_geo['gap_wide'] + ass_geo['channel_box_thickness']
            W_end = ass_geo['assembly_pitch'] - W_start

            z_coords = np.linspace(z2, z1, len(matrice_numpy))
            step_z = max(1, len(z_coords) // 15) 
            labels_z = [f"{val:.1f}" if i % step_z == 0 else "" for i, val in enumerate(z_coords)]
            
            y_coords_phys = np.linspace(W_start, W_end, len(matrice_numpy[0]))
            step_y = max(1, len(y_coords_phys) // 10)
            labels_y = [f"{val:.2f}" if i % step_y == 0 else "" for i, val in enumerate(y_coords_phys)]
            
            sns.heatmap(matrice_a_dessiner, annot=afficher_nombres, fmt=".2f", 
                        cmap="viridis", cbar_kws={'label': r'Porosité $\phi$ (-)'},
                        yticklabels=labels_z, xticklabels=labels_y)
            
            plt.yticks(rotation=0) 
            plt.xticks(rotation=45) 
            
            plt.title(f"Carte de Porosité 3D (En X={x1}cm) - {args.mesh_x.capitalize()} - Assm: {args.core_pos[0]},{args.core_pos[1]}")
            plt.xlabel("Axe Y (Largeur de l'assemblage en cm)")
            plt.ylabel("Axe Z (Altitude en cm)")

            nom_fichier = f"heatmap_mesh_x_{args.mesh_x}_assm{args.core_pos[0]}{args.core_pos[1]}_X_slice{x1}_Z{z1}_Z{z2}.png"
            chemin_complet = os.path.join(DOSSIER_OUTPUTS, nom_fichier)
            plt.savefig(chemin_complet, dpi=300, bbox_inches='tight')
            plt.close()

            print(f"Heatmap sauvegardée avec succès dans :\n -> {chemin_complet}\n")

    elif args.mesh_y:
        mat_p, mat_dh = analyser.execute_mesh_y(mesh_type_y, y1, z1, z2, n_z)
        print("--- Matrice de Porosité (phi) ---")
        print(np.round(np.array(mat_p), 4)) # Arrondi à 4 décimales pour la lisibilité
        
        print("\n--- Matrice de Diamètre Hydraulique (Dh) ---")
        print(np.round(np.array(mat_dh), 4))
        
        if args.plot:
            matrice_numpy = np.array(mat_p)
            plt.figure(figsize=(10, 8))
            matrice_a_dessiner = matrice_numpy[::-1]
            afficher_nombres = (len(mat_p[0]) <= 15 and len(mat_p) <= 15)    
            
            ass_geo = analyser.data_ref['ASSEMBLY_GEOMETRY']
            W_start = ass_geo['gap_wide'] + ass_geo['channel_box_thickness']
            W_end = ass_geo['assembly_pitch'] - W_start

            z_coords = np.linspace(z2, z1, len(matrice_numpy))
            step_z = max(1, len(z_coords) // 15) 
            labels_z = [f"{val:.1f}" if i % step_z == 0 else "" for i, val in enumerate(z_coords)]
            
            x_coords_phys = np.linspace(W_start, W_end, len(matrice_numpy[0]))
            step_x = max(1, len(x_coords_phys) // 10)
            labels_x = [f"{val:.2f}" if i % step_x == 0 else "" for i, val in enumerate(x_coords_phys)]
            
            sns.heatmap(matrice_a_dessiner, annot=afficher_nombres, fmt=".2f", 
                        cmap="viridis", cbar_kws={'label': r'Porosité $\phi$ (-)'},
                        yticklabels=labels_z, xticklabels=labels_x)
            
            plt.yticks(rotation=0) 
            plt.xticks(rotation=45) 
            
            plt.title(f"Carte de Porosité 3D (En Y={y1}cm) - {args.mesh_y.capitalize()} - Assm: {args.core_pos[0]},{args.core_pos[1]}")
            plt.xlabel("Axe X (Largeur de l'assemblage en cm)")
            plt.ylabel("Axe Z (Altitude en cm)")

            nom_fichier = f"heatmap_mesh_y_{args.mesh_y}_assm{args.core_pos[0]}{args.core_pos[1]}_Y_slice{y1}_Z{z1}_Z{z2}.png"
            chemin_complet = os.path.join(DOSSIER_OUTPUTS, nom_fichier)
            plt.savefig(chemin_complet, dpi=300, bbox_inches='tight')
            plt.close()

            print(f"Heatmap sauvegardée avec succès dans :\n -> {chemin_complet}\n")

    elif args.cv_z or args.rod_z or args.water_z:
        if section_type[0] == 'cv':
            x1, y1, x2, y2 = section_type[1]
            phi = analyser.get_porosity_z_cv(x1, y1, x2, y2, z1, z2)
            a_cool = analyser.get_a_cool_z_cv(x1, y1, x2, y2, z1, z2)
            dh = analyser.get_dh_z_cv(x1, y1, x2, y2, z1, z2)
            ph = analyser.get_ph_cv(x1, y1, x2, y2, z1, z2)
            p_box = analyser.get_pbox_cv(x1, y1, x2, y2, z1, z2)
            p_wr = analyser.get_pwr_cv(x1, y1, x2, y2, z1, z2)
        elif section_type[0] == 'rod':
            i, j = section_type[1]
            phi = analyser.get_porosity_z_rod(i, j, z1, z2)
            a_cool = analyser.get_a_cool_z_rod(i, j, z1, z2)
            dh = analyser.get_dh_z_rod(i, j, z1, z2)
            ph = analyser.get_ph_rod(i, j, z1, z2)
            p_box = analyser.get_pbox_rod(i, j, z1, z2)
            p_wr = analyser.get_pwr_rod(i, j, z1, z2)
        elif section_type[0] == 'water':
            i, j = section_type[1]
            phi = analyser.get_porosity_z_canal(i, j, z1, z2)
            a_cool = analyser.get_a_cool_z_canal(i, j, z1, z2)
            dh = analyser.get_dh_z_canal(i, j, z1, z2)
            ph = analyser.get_ph_canal(i, j, z1, z2)
            p_box = analyser.get_pbox_canal(i, j, z1, z2)
            p_wr = analyser.get_pwr_canal(i, j, z1, z2)
        print(f"Porosité (phi) : {phi}")
        print(f"Surface de passage fluide (A_cool) : {a_cool} cm²")
        print(f"Dh (cm)        : {dh}")
        print(f"Ph (cm)        : {ph}\n")
        print(f"P_box (cm) : {p_box}")
        print(f"P_wr (cm) : {p_wr}")
    
    elif args.cv_x or args.rod_x or args.water_x:
        if section_type[0] == 'cv':
            y1, z1, y2, z2 = section_type[1]
            phi = analyser.get_porosity_x_cv(x1, y1, y2, z1, z2)
            dh = analyser.get_dh_x_cv(x1, y1, y2, z1, z2)
        elif section_type[0] == 'rod':
            i, j = section_type[1]
            phi = analyser.get_porosity_x_rod(i, x1, z1, z2)
            dh = analyser.get_dh_x_rod(i, x1, z1, z2)
        elif section_type[0] == 'water':
            i, j = section_type[1]
            phi = analyser.get_porosity_x_canal(i, x1, z1, z2)
            dh = analyser.get_dh_x_canal(i, x1, z1, z2)
        print(f"Porosité (phi) : {phi}")
        print(f"Dh (cm)        : {dh}")

    elif args.cv_y or args.rod_y or args.water_y:
        if section_type[0] == 'cv':
            x1, z1, x2, z2 = section_type[1]
            phi = analyser.get_porosity_y_cv(y1, x1, x2, z1, z2)
            dh = analyser.get_dh_y_cv(y1, x1, x2, z1, z2)
        elif section_type[0] == 'rod':
            i, j = section_type[1]
            phi = analyser.get_porosity_y_rod(j, y1, z1, z2)
            dh = analyser.get_dh_y_rod(j, y1, z1, z2)
        elif section_type[0] == 'water':
            i, j = section_type[1]
            phi = analyser.get_porosity_y_canal(j, y1, z1, z2)
            dh = analyser.get_dh_y_canal(j, y1, z1, z2)
        print(f"Porosité (phi) : {phi}")
        print(f"Dh (cm)        : {dh}")
        
        

    
    

