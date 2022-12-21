
<h4>This utility is used together with buoyantPimpleFoam to calculate the wall Nusselt number.</h4>

<h3>Usage:</h5>

<ol>

<li>controlDict:</li>
<pre>
    wallNu
    {
        type            myWallNu;
        libs            ( "libNu_rho_v2.so" );
        writeControl    writeTime;
        Tref            300; //reference temperature;
        Lref            0.005; //reference length; 
    }
</pre>

</ol>
