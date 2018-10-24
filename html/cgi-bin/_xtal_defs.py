"""
defs for xtal.cgi
"""

ERR_HTML = """\
Content-Type: text/html\n
<html>
<body>
<p>
%s
</p>
</body>
</html>
"""

###########

START_HEAD = """<!DOCTYPE html>
<head>
<title>Xtal Output</title>
<meta charset="utf-8">
"""

JMOL_1_SCRIPT = """
<!--  Jmol script for structure.xyz -->

<script type="text/javascript" src="../../jmol/JSmol.min.js"></script>

<script type="text/javascript"> 

$(document).ready(function() {

script = 'load structure.xyz;'
+ 'boundbox; axes 1; set perspectiveDepth OFF; '
+ 'connect (all) (all) Delete; connect 2.2 (all) (all) Create; '
+ 'wireframe 0.1; spacefill 20%;'
+ 'moveto 0 AXIS x; '

Info = {
    debug: false,
    color: "0xC0C0C0",
    height: 600,
    width: 800,
    j2sPath: "../../jmol/j2s",
    jarPath: "../../jmol/java",
    jarFile: "JmolAppletSigned0.jar",
    isSigned: false,
    serverURL: "../../jmol/php/jsmol.php",
    disableJ2SLoadMonitor: true,
    disableInitialConsole: true,
    addSelectionOptions: false,
    use: "HTML5",
    readyFunction: null,
    script: script
}

$("#struct_xyz").html(Jmol.getAppletHtml("jmolApplet0",Info))

});

</script>
"""

JMOL_2_SCRIPT = """
<!--  Jmol script for bulk.cif showing HKL plane -->

<script type="text/javascript"> 

$$(document).ready(function() {

script = 'load bulk.cif {1 1 1}; '
+ 'axes 3;'
+ 'connect (all) (all) Delete; connect 2.2 (all) (all) Create;'
+ 'wireframe 0.1; spacefill 20%;'
+ 'moveto 0 AXIS a4; '
+ 'isosurface hkl {$H  $K  $L}  translucent'

Info = {
    debug: false,
    color: "0xC0C0C0",
    height: 600,
    width: 800,
    j2sPath: "../../jmol/j2s",
    jarPath: "../../jmol/java",
    jarFile: "JmolAppletSigned0.jar",
    isSigned: false,
    serverURL: "../../jmol/php/jsmol.php",
    disableJ2SLoadMonitor: true,
    disableInitialConsole: true,
    addSelectionOptions: false,
    use: "HTML5",
    readyFunction: null,
    script: script
}

$$("#bulk_hkl").html(Jmol.getAppletHtml("jmolApplet1",Info))

});

</script>
"""

END_HEAD = """</head>
"""

START_HTML = """
<!-----------------  start html  --------------------->

<body>

<h1>Xtal Stucture Output:  ${FNAME}  </h1>

<table border="1">
"""

XTAL_1_OUT = """
  <!-- xtal output-->
  <tr>
    <td>
        <h3> Structure information: <a href="structure.xyz" target="_blank">  Download structure.xyz</a></h3>  
        <textarea name="bulk_xtal" cols="130" rows="25">  ${XTAL} </textarea>
        </br></br>
    </td>
  </tr>
"""

XTAL_2_OUT = """
  <!-- xtal output-->
  <tr>
    <td>
        <h3> Structure information: <a href="structure.xyz" target="_blank">  Download structure.xyz</a> 
             <a href="structure.frac" target="_blank">  Download structure.frac</a></h3>  
        <textarea name="bulk_xtal" cols="130" rows="25">  ${XTAL} </textarea>
        </br></br>
    </td>
  </tr>
"""

COORD_OUT = """
  <!-- Coord calcs-->
  <tr>
    <td>
        <h3>Coord calculations: <a href="coord.out" target="_blank">  Download coord.out</a></h3>
        <textarea name="coord" cols="130" rows="25">  ${COORD} </textarea>
    </td>
  </tr>
"""

STRUCT_1_OUT = """
  <!-- Structure-->
  <tr>
    <td>
        <h3>Jmol: structure.xyz </h3>
        <span id=struct_xyz></span>
    </td>
  </tr>
"""

STRUCT_2_OUT = """
  <!-- Structure-->
  <tr>
    <td>
        <h3>Jmol: bulk xtal showing {$H $K $L}</h3>
        <span id=bulk_hkl></span>
    </td>
  </tr>
"""


END = """
</table></body></html>
"""

if __name__ == "__main__":
    html = START_HEAD
    html = html + JMOL_SCRIPT
    html = html + END_HEAD
    html = html + START_HTML
    html = html + XTAL_OUT
    html = html + COORD_OUT
    html = html + STRUCT_1_OUT
    html = html + END
    f = open('test.html','w')
    f.write(html)
    f.close()
    print(html)





