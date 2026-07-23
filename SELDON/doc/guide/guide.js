
function DisplayAccor()
{
    var acc = document.getElementsByClassName("accordeon");
    var i;

    var chaine = window.localStorage.getItem("SeldonMenu");
    if (chaine != null)
    {
        var mots = chaine.split(" ");
        if (mots.length >= acc.length)
        {
            for (i = 0; i < acc.length; i++)
            {
                if (mots[i] == "1")
                {
                    acc[i].classList.add("activeAccor");
                    var nom_panel = acc[i].getAttribute("name");
                    var panel = document.getElementById(nom_panel);
                    panel.style.display = "block";
                }
            }
        }
    }
    
    for (i = 0; i < acc.length; i++) {
        acc[i].addEventListener("click", function() {
            /* Toggle between adding and removing the "active" class,
               to highlight the button that controls the panel */
            this.classList.toggle("activeAccor");
            var j; var chaine = "";
            for (j = 0; j < acc.length; j++)
            {
                if (acc[j].classList.contains("activeAccor"))
                    chaine += "1 ";
                else
                    chaine += "0 ";
            }
            
            window.localStorage.setItem("SeldonMenu", chaine);
            
            /* Toggle between hiding and showing the active panel */
            var nom_panel = this.getAttribute("name");
            var panel = document.getElementById(nom_panel);
            if (panel.style.display === "block") {
                panel.style.display = "none";
            } else {
                panel.style.display = "block";
            }
        });
    }
}
