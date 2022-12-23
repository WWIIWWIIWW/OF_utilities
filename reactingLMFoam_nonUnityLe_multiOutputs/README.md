
<h4>This solver is based on <a href="https://github.com/WWIIWWIIWW/OF_utilities/tree/main/reactingLMFoam_nonUnityLe">reactingLMFOAM_nonUnityLe</a>.</h4>

<h3>Instructions:</h5>

<ol>

<li>LESResIndex.H:</li>
<pre>
    Calculate and update many variables which is why this solver is called multiOutputs.
    (we can have some written in "writeObj1", but not all.)
    Calculated kRes and kTOT, the two are used for calculating LESIndex later.
</pre>

<li>updateVariable.H:</li>
<pre>
    Variables updated in this folder was claimed in "createFieldRefs.H".
    AUTO_WRITE Enabled.
    Eg:
    <hr>
    volScalarField Cp_
    (
        IOobject
        (
            "Cp_",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            <b>IOobject::AUTO_WRITE</b>
        ),
        thermo.Cp()
    );
</pre>
</ol>
