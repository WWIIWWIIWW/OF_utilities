IOdictionary IOsnapEdgeDict
(
    IOobject
    (
        "snapEdgeDict",
         runTime.constant(),
         runTime,
         IOobject::MUST_READ,
         IOobject::NO_WRITE
    )
);

dictionary snapEdgeDict(IOsnapEdgeDict.subDict("snapEdgeDict"));
