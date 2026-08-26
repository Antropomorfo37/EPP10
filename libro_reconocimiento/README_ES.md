# Libro de reconocimiento

Aplicación **local y autocontenida** para digitalizar las notas de reconocimiento
que los alumnos escriben al final de las clases, organizarlas en capítulos por
grupo y componer con ellas un libro imprimible.

- **Todo queda en tu dispositivo.** La app es un único archivo (`index.html`)
  que funciona sin internet; las fotos y los datos se guardan en el propio
  navegador (IndexedDB). Nada se sube a ningún servidor y nada se guarda en
  este repositorio.
- **Privacidad por defecto.** En la app y en el libro no aparece ningún nombre.
  Solo existe una versión especial del libro **con nombres**, pensada
  únicamente para entregarla al jefe, que hay que activar de forma explícita
  en cada exportación.

## Cómo abrirla

- **En el ordenador (recomendado):** haz doble clic en
  `libro_reconocimiento/index.html` (Chrome, Edge o Firefox).
- **Alternativa 100 % local por navegador estricto:** desde la carpeta del
  repositorio ejecuta `python3 -m http.server 8000` y abre
  `http://localhost:8000/libro_reconocimiento/`.

*(Esta guía se completará con el flujo de trabajo detallado — captura por
tandas, tapado de nombres, anexos, exportación a PDF y copias de seguridad —
a medida que se añadan las funciones.)*
