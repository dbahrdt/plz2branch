## Eine neue Fapra-Aufgabe aus der realen Welt

# Input File Format

PLZ is just a hint
```
<NAME-1> <LATITUDE-1> <LONGITUDE-1> <PLZ1> <PLZ2> <PLZ3>
<NAME-2> <LATITUDE-2> <LONGITUDE-2> <PLZ1> <PLZ2> <PLZ3>
```

## Distanz
Jeder PLZ wird die jeweils naheste Filiale zugeordnet.
Nahe bedeutet hier im Regelfall die mittlere Distanz.
Nahe kann hierbei noch gewichtet werden:
 - Mit der Mitarbeiteranzahl der Filiale
 - Der Bevölkerungsgröße der PLZ
 
## Greedy pick
PLZ verteilen in Runden:
In jeder Runde wird für jede Filiale gewählt, die die Kosten der Filiale am wenigsten erhöht.
Diejenige Filiale mit den zukünftigen geringsten Gesamtkosten bekommt den Zuschlag.

Vorgehen bei bereits zugewiesen PLZ
1. Überspringen
2. Wegnehmen

Falls hierdurch eine PLZ einer anderen Filiale weggenommen wird, dann werden diese Kosten der entsprechenden Filiale wieder entfernt.
Um ein Oszilieren zu vermeiden steigen die Kosten für das Wählen einer PLZ mit jedem Wählen dieser PLZ an.

## Evolutionär
DNA ist der Vektor V der Länge |P| mit Einträgen aus F.
Mutation:
 - Beliebigen Eintrag durch ein random Element aus F ersetzen
Rekombination:
 - V_1, V_2 linear durchlaufen und jeweils den Eintrag nehmen, der die geringeren Kosten hat
Selektion:
Population nach Kosten sortieren und aus jedem Quantil exponentiell absteigend viele wählen.
D.h. von den "guten" gibt es viel mehr als von den "schlechten".
Mit der Zeit sollten die Kosten der "schlechten" kleiner werden.

## Optimal Problem Formulation

$F = \left\{ f_0 \ldots f_m \right\}$ Menge an Filialen
$P = \left\{ p_0 \ldots p_n \right\}$ Menge an Postleitzahlen
$c: F \times P \to \R, (f, p) \mapsto Distanz(f, p) \cdot \frac{Bevölkerung(p)}{Mitarbeiter(f)}$ Kostenfunktion
Gesucht ist Menge $E \subset F \times P$ mit $cost = \Sigma_{(u,v) \in E c(u,v)}$ minimal und $P = \{ p \in P | (f,p) \in E\}$
### Min-Cost-Max-Flow
Geht. Bild folgt.

### Evolutionärer Algorithmus
Geht einfach.
