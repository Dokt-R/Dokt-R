/*
Ashely(P)
Christeen(P)
Jane(A)
Jenny(D)
Julia(A)
Ketty(P)
Maria(A)
Meera(S)
Priya(S)
Samantha(D)
There are a total of 2 doctors.
There are a total of 2 singers.
There are a total of 3 actors.
There are a total of 3 professors.
*/
SELECT CONCAT(name,'(',left(occupation,1),')') FROM occupations ORDER BY name ASC;
SELECT CONCAT('There are a total of ',count(occupation),' ',lower(occupation),'s.') FROM occupations GROUP BY occupation ORDER BY count(occupation), occupation ASC