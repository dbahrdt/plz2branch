# Input File Format
```
<NAME-1>
<LATITUDE-1>
<LONGITUDE-1>
<NAME-2>
<LATITUDE-2>
<LONGITUDE-2>
```

## What we really want?

```C++

Sei P={p_0, ... p_n} eine Menge von Postleitzahlen.
F={f_0, ... f_m} eine Menge von Filialen.
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

```
