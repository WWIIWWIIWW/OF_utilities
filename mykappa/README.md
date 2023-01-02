
<h4>This utility is used together with reactingLMFoam_nonUnityLe_multiOutputs to calculate the near-wall kappa (1 / (1+Da)).</h4>

<h3>Usage:</h5>

<ol>

<li>controlDict:</li>
<pre>
    kappa
    {
        type            mykappa;
        libs            ( "libkappa.so" );
        //executeControl  timeStep;
        //executeInterval 1;
        writeControl    timeStep;//writeTime;
        writeInterval   1;
    }
</pre>

</ol>
