import os, sys

# pour inserer le menu dans le fichier html
def insert_menu(fichier_html, fichier_menu):
    f = open(fichier_html, "r")
    lignes_html = f.readlines()
    f.close()

    f = open(fichier_menu, "r")
    lignes_menu = f.readlines()
    f.close()

    pos_pagedoc = 0; pos_contents = 0
    for i in range(len(lignes_html)):
        if (lignes_html[i].startswith('<div class="PageDoc">')):
            pos_pagedoc = i
            
        if (lignes_html[i].startswith('<div class="contents"')):
            pos_contents = i

    # on enleve la ligne avant PageDoc avec </div><!-- top -->
    for i in range(pos_pagedoc-1, pos_contents):
        lignes_html.pop(pos_pagedoc-1)

    # puis on insere le menu
    for i in range(len(lignes_menu)):
        lignes_html.insert(pos_pagedoc-1+i, lignes_menu[i])

    for i in range(len(lignes_html)):
        if (lignes_html[i].startswith('</div> <!-- doxygen -->')):
            lignes_html.pop(i)
            break
        
    # et on ecrase le fichier
    f = open(fichier_html, "w")
    for i in range(len(lignes_html)):
        f.write(lignes_html[i])
    
    f.close()

# pour enlever certains </div> parasites
def remove_div(fichier_html):
    f = open(fichier_html, "r")
    lignes_html = f.readlines()
    f.close()

    for i in range(len(lignes_html)):
        if (lignes_html[i].startswith('</div><!-- top -->')):
            lignes_html.pop(i)
        
        if (lignes_html[i].startswith('</div> <!-- doxygen -->')):
            lignes_html.pop(i)
            break
    
    # et on ecrase le fichier
    f = open(fichier_html, "w")
    for i in range(len(lignes_html)):
        f.write(lignes_html[i])
    
    f.close()

# pour bien mettre prettyprint
def clean_prettyprint(fichier_html):
    f = open(fichier_html, "r")
    lignes_html = f.readlines()
    f.close()

    chaine = '<pre class="prettyprint"><pre class="fragment">'
    chaineF = '</pre> </pre>';
    lc = len(chaine); lcF = len(chaineF);
    modif = False; i = 0;
    while (i < len(lignes_html)):
        if (lignes_html[i].find(chaine) >= 0):
            pos = lignes_html[i].find(chaine)
            p = lignes_html[i];
            modif = True
            lignes_html[i] = p[0:pos] + '<pre class="prettyprint">' + p[pos+lc:]
            while (lignes_html[i].find(chaineF) < 0):
                i += 1;

            pos = lignes_html[i].find(chaineF)
            p = lignes_html[i];
            lignes_html[i] = p[0:pos] + '</pre>' + p[pos+lcF:]
        
        i += 1;
    
    # et on ecrase le fichier si modif
    if (modif):
        f = open(fichier_html, "w")
        for i in range(len(lignes_html)):
            f.write(lignes_html[i])
    
        f.close()
    
liste_fichiers = sorted(os.listdir("."))
for f in liste_fichiers:
    if (f.endswith('.dox')):
        racine = f[0:f.find('.dox')]
        nom_menu = "menu_" + racine
        nom_html = "../html/" + racine + ".html"
        if (racine == "guide"):
            nom_html = "../html/index.html"
        if (nom_menu in liste_fichiers):
            insert_menu(nom_html, nom_menu)
        else:
            remove_div(nom_html)

        clean_prettyprint(nom_html)



