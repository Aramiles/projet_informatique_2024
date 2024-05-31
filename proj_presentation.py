import streamlit as st
from projet_info import molecule
from rdkit import Chem
from rdkit.Chem import Draw
import os
import pandas as pd
import time
import re
from selenium import webdriver
from selenium.webdriver.firefox.service import Service
from selenium.webdriver.common.by import By
from selenium.webdriver.common.keys import Keys
from selenium.webdriver.firefox.options import Options

def wiki_search(name):
    firefox_driver_path = "geckodriver-v0.33.0-win64/geckodriver.exe"
    firefox_options = Options()
    firefox_options.add_argument('--headless')  # Activer le mode headless
    
    service = Service(executable_path=firefox_driver_path)
    browser = webdriver.Firefox(service=service, options=firefox_options)
    browser.get("https://en.wikipedia.org/wiki/Main_Page")
    
    # Localiser la zone de recherche et saisir le texte
    search_box = browser.find_element(By.ID, 'searchInput')
    search_box.send_keys(name)
    search_box.submit()
    
    time.sleep(1)
    
    paragraphs = browser.find_elements(By.TAG_NAME, 'p')
    first_paragraph_with_text = None
    for paragraph in paragraphs:
        if paragraph.text.strip() != "":
            first_paragraph_with_text = paragraph
            break
    
    clean_text = re.sub(r'\[\d+\]', '', first_paragraph_with_text.text)
    st.write(clean_text)
    st.write(f"from the english wikipedia page about {name}")
    
    browser.quit()


    
            #title

st.title("project of molecule's analysis")

st.subheader("1)   **let's charge datas**")

            #charge datas

user_input = st.text_input("write the name of your folder on this case")
try:
    new_mol = molecule(user_input)
    if new_mol.load_verify_st() == True :
        mole= molecule(user_input)
        # Afficher la saisie de l'utilisateur
        st.write("You want to analyse the molecule: ", user_input[:-4])
        number_of_atoms = len(mole.get_atoms)

        button_clicked_mol = st.button("small description")
        if button_clicked_mol:
            wiki_search(user_input[:-4])


        
        
        col1, col2 = st.columns(2)

           #show molecules

        
        with col1:
            button_clicked_mol = st.button("show molecule")
            if button_clicked_mol:
                img_mole = Draw.MolToImage(Chem.MolFromMolFile(os.path.join(os.path.dirname(os.path.abspath(__file__))+"\MOL_files", user_input)))
                st.image(img_mole, caption="Molecular structure", width=300)     
                
        with col2:
                st.table(pd.DataFrame(mole.get_infos.values(), index = mole.get_infos.keys(), columns = ["datas"]))
        
        
        st.markdown("""<style>.centered-line {width: 90%;margin: 0 auto;border: 1px solid #ddd}</style>""", unsafe_allow_html=True)
        st.markdown("<hr class='centered-line'>", unsafe_allow_html=True)

        
        col1, col2 = st.columns(2)


        
        # show .mol file

        
        with col1:
            st.subheader("about atoms")
            button_clicked_atom = st.button("show atom's bloc")
        
        if button_clicked_atom:
            st.table(mole.get_atoms)
        

        
        with col2:
            st.subheader("about bonds")
            button_clicked_bonds = st.button("show bond's bloc")
        
            # Afficher le tableau si le bouton est cliqué
        if button_clicked_bonds:
            st.table(mole.get_bonds)
        
        st.markdown("""<style>.centered-line {width: 90%;margin: 0 auto;border: 1px solid #ddd}</style>""", unsafe_allow_html=True)
        st.markdown("<hr class='centered-line'>", unsafe_allow_html=True)
        
        st.subheader("2)   **Some functionalities**")


            #project's functions

        
        col1, col2 = st.columns(2)



            #count elements
        
        with col1:
            st.subheader("A- count elements")
            st.subheader("")
            element = st.text_input("Please enter the type of atom by its letter")
            if mole.count_elements(element) != 0:
                st.success(f"The number of {element} atom elements in the molecule is {mole.count_elements(element)}.") 
            elif element == "":
                st.warning("No element provided")
            else:
                st.warning(f"The atom element {element} does not exist in the molecule.")



                 #2D distance
            
            st.subheader("C- Topological distance between two atoms of the molecule")
            couple = st.text_input(f"Please enter the identifier of atom 1 and 2 with a space between (between 1 and {number_of_atoms}) ")
            
            try:
                
                atom1_top, atom2_top = map(int, couple.split())
                
                if atom1_top == atom2_top:
                    st.success(f"The distance between the atoms {atom1_top} and {atom2_top} in the molecule is 0.")
                elif (1 <= atom1_top <= number_of_atoms) and (1 <= atom2_top <= number_of_atoms):
                    distance = mole.distance_topo_2D(atom1_top, atom2_top)
                    st.success(f"The distance between the atoms {atom1_top} and {atom2_top} in the molecule is {distance}.")
                else:
                    st.warning(f"One or both of the atoms {atom1_top} and {atom2_top} do not exist in the molecule.")
            except ValueError:
                st.warning("Please enter two valid atom identifiers separated by a space.")
            except Exception as e:
                st.warning(f"An error occurred: {e}")


                #3D distance

        
        with col2:
            
            # 3D distance
            st.subheader("B- 3D Distance Between Two Atoms of the Molecule")
            couple3D = st.text_input(f"Please enter the identifier of atom 1 and 2 with a space between (between 1 and {number_of_atoms})")
            
            try:
                atom1_3D, atom2_3D = map(int, couple3D.split())
                if atom1_3D == atom2_3D:
                    st.success(f"The distance between the atoms {atom1_3D} and {atom2_3D} in the molecule is 0.")
                elif 1 <= atom1_3D <= number_of_atoms and 1 <= atom2_3D <= number_of_atoms:
                    distance3D = mole.distance_3D(atom1_3D, atom2_3D)
                    st.success(f"The distance between the atoms {atom1_3D} and {atom2_3D} in the molecule is {distance3D}.")
                else:
                    st.warning(f"One or both of the atoms {atom1_3D} and {atom2_3D} do not exist in the molecule.")
            except ValueError:
                st.warning("Please enter two valid atom identifiers separated by a space.")
            except Exception as e:
                st.warning(f"An error occurred: {e}")
        
            # Neighbours
            st.subheader("D- Atom's Neighbours")
            ato = st.text_input(f"Please enter the identifier of an atom (between 1 and {number_of_atoms})")
            
            if ato.strip() == "":
                st.warning("No atom provided")
            else:
                try:
                    ato = int(ato)
                    if 1 <= ato <= number_of_atoms:
                        neighbours = mole.atom_neighbours(ato)
                        if neighbours:
                            st.success(f"The number of neighbours for atom {ato} is {len(neighbours)}.")
                        else:
                            st.warning(f"The atom {ato} has no neighbours.")
                    else:
                        st.warning(f"The atom {ato} does not exist in the molecule.")
                except ValueError:
                    st.warning("Please enter a valid integer for the atom identifier.")
                except Exception as e:
                    st.warning(f"An error occurred: {e}")


                #ring
                    
        st.subheader("E- Ring verification")
        st.subheader("")
        ring_chain = st.text_input(f"Please enter the identifiers of atoms you want to check if they form a ring. Put an espace between each indentifier (identifiers should be between 1 and {number_of_atoms}) ")
        mole.ring_finding_st(ring_chain)

    st.subheader("3)   **Conclusion**")
    st.write(f"It is possible and simple to create and add new features without significantly modifying the source code. Students in the next years will be able to make this page more complete. We hope we were able to help you with your file: {user_input}")
        
except:
    st.write("File not found or no file provided")


