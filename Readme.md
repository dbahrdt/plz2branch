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
$c: F \times P \to \R, (f, p) \mapsto dist(f, p)*people(p)/workers(f)$ Kosten, falls Filiale $f_i$ die Postleitzahl $p_j$ bedient mit $dist$ einer Funktion, die die Distanz von der Filiale zur PLZ angibt.
Gesucht ist eine Überdeckung von $P$ mit Teilmengen $Z_i \subseteq P$ mit $\forall i,j : i \neq j \Rightarrow Z_i \cap Z_j = \emptyset$ und $cost = \Sigma_i Z_i = \Sigma_{u,v \in F \times P} c(f, p)$
### Min-Cost-Max-Flow
$S_F = \{(u, v) : u = s \wedge v \in F \}$
E=S

```

### ILP

Sei $P={p_0, ... p_n}$ eine Menge von Postleitzahlen.
$F=\{f_0, ... f_m\}$ eine Menge von Filialen.
m: F -> \N eine Funktion, die die Mitarbeiternanzahl einer Filiale widergibt
b: P -> \R eine Funktion, die die Bevölkerungsdichte einer Postleitzahl angibt
d: P x F -> \R eine Funtion, die die mittlere Distanz einer Postleitzahl zu einer Filiale engibt.
Dann seien c_ij=d(p_i, f_j)*b(p_i)/m(f_j) die Kosten der Filiale f_j für die PLZ p_i.
Desweiteren sei x_ij \in {0, 1} eine Indikator-Variable, die angibt, ob PLZ p_i von Filiale f_j bedient wird.
Weiterhin sein l_j \in \R eine Konstante, die angibt, wieviele PLZ eine Filiale mindestens bedienen muss.
Dabei gilt \sum_j l_j >= |P|.
Weiterhin sei l_j <= u_j \in \R eine Konstante, die angibt, wieviele PLZ eine Filiale maximal bedienen muss.

Dann wollen wir folgendes Optimierungsproblem lösen:

min \sum_{i,j} c_ij * x_ij
s.t.

x_ij \in {0, 1}

Für jede PLZ p_i exakt eine Filiale f_j
\sum_j x_ij = 1

Für jede Filiale f_j mindestens l_j viele PLZ
    l_j <= \sum_i x_ij
<=> \sum_i -x_ij <= -l_j

Für jede Filiale f_j maximal u_j viele PLZ
\sum_i x_ij <= u_j
