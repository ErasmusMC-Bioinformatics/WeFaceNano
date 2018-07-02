function openResult(evt, rname) {
    // Declare all variables
    var i, tabcontent, tablinks;
    // Get all elements with class="tabcontent" and hide them
    tabcontent = document.getElementsByClassName("tabcontent");
    for (i = 0; i < tabcontent.length; i++) {
        tabcontent[i].style.display = "none";
    }
    // Get all elements with class="tablinks" and remove the class "active"
    tablinks = document.getElementsByClassName("tablinks");
    for (i = 0; i < tablinks.length; i++) {
        tablinks[i].className = tablinks[i].className.replace(" active", "");
    }
    // Show the current tab, and add an "active" class to the button that opened the tab
    document.getElementById(rname).style.display = "block";
    evt.currentTarget.className += " active";
};
function openTab(evt, rname) {
    // Declare all variables
    var i, tabdiv, tabbutton;
    // Get all elements with class="tabcontent" and hide them
    tabdiv = document.getElementsByClassName("tabdiv");
    for (i = 0; i < tabdiv.length; i++) {
        tabdiv[i].style.display = "none";
    }
    // Get all elements with class="tablinks" and remove the class "active"
    tabbutton = document.getElementsByClassName("tabbutton");
    for (i = 0; i < tabbutton.length; i++) {
        tabbutton[i].className = tabbutton[i].className.replace(" active", "");
    }
    // Show the current tab, and add an "active" class to the button that opened the tab
    document.getElementById(rname).style.display = "block";
    evt.currentTarget.className += " active";
};
