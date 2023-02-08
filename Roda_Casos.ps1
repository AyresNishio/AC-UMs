$arquivos =  Get-ChildItem ".\Casos" -Filter *.txt
Write-Output $arquivos
g++ .\CritFindMU.cpp -o critfind
foreach($caso in $arquivos)
{
    Copy-Item $caso.FullName -Destination .\Caso.txt
    .\critfind.exe
    $saida = $caso.BaseName+".csv"
    Copy-Item .\Saida.csv -Destination .\Casos\$saida
}