# vim: set fileencoding=utf-8 :

def generate_html_table(mibig_homology_file):
    outfile_name  = mibig_homology_file.split('.txt')[0] + '.html'
    outfile = open(outfile_name,"w")
    outfile.write("""

<!DOCTYPE html PUBLIC "-//W3C//DTD XHTML 1.0 Transitional//EN" "http://www.w3.org/TR/xhtml1/DTD/xhtml1-transitional.dtd">
<html xmlns="http://www.w3.org/1999/xhtml">

<head>
<link rel="stylesheet" type="text/css" href="../../css/style.css">
<!-- js section -->
<!-- jQuery  -->
<script type="text/javascript"
        src="../../js/jquery-1.10.2.min.js"></script>

<script type="text/javascript" src="../../js/jquery.tablesorter.min.js"></script>

<script type="text/javascript" src="../../js/jquery.tmpl-1.0.0pre.js"></script>

  <!-- jQuery UI Support -->
<script type="text/javascript"
        src="../../js/jquery-ui-1.10.4.custom.min.js"></script>

 <!-- Alpaca -->
<script type="text/javascript" src="../../js/alpaca-full-mibig.js"></script>

<script type="text/javascript" >
var counter = 1;
function addInput(divName){
     var newdiv = document.createElement('div');
     newdiv.innerHTML = "<b>Locus " + (counter + 1) + "</b><br> GenBank accession: <input type='text' required name='myInputs[]'> Start nt: <input type='number' name='myStartCoords[]' placeholder='N/A'> End nt: <input type='number' name='myEndCoords[]' placeholder='N/A'> <br><br>";
     document.getElementById(divName).appendChild(newdiv);
     counter++;
}
</script>

<script>
$(document).ready(function()
    {
        $("#myTable").tablesorter();
    }
);

</script>
<style>
table {
    border-width: thin;
    border-spacing: 2px;
    border-style: none;
    border-color: black;
}
</style>
</head>
    <table id="myTable" class="tablesorter" border=1>
    <thead>
    <tr>
        <th>MIBiG Protein</th>
        <th>Description</th>
        <th>MIBiG Cluster</th>
        <th>MiBiG Product</th>
        <th>% ID</th>
        <th>% Coverage</th>
        <th>BLAST Score</th>
        <th>E-value</th>
    </tr>
    </thead>
    <tbody> """)
    for line in open(mibig_homology_file):
        if line[0] == '#':
            pass
        else:
            line_parse = line.split('\t')
            line_parse = map(lambda x:x.strip(),line_parse)
            gene_id,gene_descr,mibig_cluster,mibig_product,percent_id,percent_cvg,blast_score,e_value = line_parse
            mibig_id = mibig_cluster.split('_c1')[0]
            outfile.write("""
            <tr>
                <td>%s</td>
                <td>%s</td>
                <td><a href="http://mibig.secondarymetabolites.org/repository/%s/index.html#cluster-1">%s</a></td>
                <td>%s</td>
                <td>%s</td>
                <td>%s</td>
                <td>%s</td>
                <td>%s</td>
            </tr>
            """ % (gene_id, gene_descr, mibig_id, mibig_id,mibig_product, percent_id, percent_cvg,blast_score,e_value))

    outfile.write("""
    </tbody>
    </table>
    """)