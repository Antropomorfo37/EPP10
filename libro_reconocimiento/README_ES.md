# Libro de reconocimiento

Aplicación **local y autocontenida** para digitalizar las notas de reconocimiento
que los alumnos escriben al final de las clases, organizarlas en capítulos por
grupo y componer con ellas un **libro imprimible** — en dos versiones: una
anónima (la normal) y otra con nombres, solo para entregar al jefe.

Toda la aplicación es un único archivo, `index.html`. No necesita internet, ni
servidor, ni instalación. Las fotos y los datos se guardan en el propio
navegador (IndexedDB) y **nunca salen del dispositivo**.

---

## 1. Privacidad primero

- **Por defecto no aparece ningún nombre**, ni en la aplicación ni en el libro.
  El nombre de cada alumno se guarda como dato privado (dentro de un desplegable
  cerrado en cada nota) y solo se usa en la versión con nombres y en tus copias
  de seguridad.
- Como las notas manuscritas suelen ir **firmadas**, la aplicación incluye una
  herramienta para **tapar** con rectángulos negros los nombres visibles en la
  foto. Se conserva la foto original y se genera una copia tapada: la versión
  anónima del libro usa la copia tapada; la versión para el jefe usa la original.
- La **lista de clase** y las **valoraciones** (anexos) solo se incluyen en la
  versión con nombres.
- La versión con nombres exige marcar una casilla de confirmación en cada
  exportación, lleva la marca «CONFIDENCIAL» en cada página y **nunca queda
  como estado guardado**: al recargar, el libro siempre vuelve a la versión
  anónima.

> Estás manejando datos personales de menores. Guarda las copias de seguridad
> (que sí contienen nombres) en lugar seguro y entrega la versión con nombres
> únicamente a quien corresponda.

---

## 2. Cómo abrir la aplicación

### En el ordenador (recomendado)
Haz doble clic en `libro_reconocimiento/index.html`. Se abre en tu navegador
(Chrome, Edge o Firefox). Funciona sin conexión.

### Si tu navegador es estricto con `file://`
Algunos navegadores (sobre todo **Safari**) restringen el almacenamiento cuando
la página se abre como archivo. En ese caso, sírvela en local — sigue siendo
100 % local, nada sale a internet:

```bash
# desde la carpeta del repositorio
python3 -m http.server 8000
```

y abre `http://localhost:8000/libro_reconocimiento/`.
(En Windows/macOS con Python instalado, el comando es el mismo.)

### En el móvil
- **Android:** puedes usar la aplicación directamente y el botón «Hacer foto»
  abrirá la cámara.
- **iPhone:** Safari no abre bien archivos `file://` con JavaScript. Lo más
  cómodo es **hacer las fotos con la cámara normal del iPhone**, pasarlas al
  ordenador (AirDrop, iCloud o cable) y allí importarlas con «Elegir fotos».
  Si prefieres trabajar en el propio iPhone, usa la alternativa `localhost` de
  arriba.

---

## 3. Flujo de trabajo recomendado

### Digitalizar tu colección histórica (desde ~2019)
1. Crea un **grupo** por cada clase o generación (p. ej. «Generación 2019»).
   Ordena los grupos cronológicamente con las flechas ↑ ↓: serán los capítulos.
2. Abre el grupo. En la barra de captura, elige la **precisión de la fecha**:
   para tandas antiguas, «Solo año» (2019) o «Sin fecha» — el capítulo ya aporta
   el contexto temporal.
3. Pulsa «**Hacer foto**» para fotografiar varias notas seguidas, o
   «**Elegir fotos**» para importar de golpe las que ya tengas escaneadas.
   Cada foto se convierte al instante en una nota **pendiente**.
4. Cuando termines la tanda, ve completando las pendientes con calma (filtro
   «Pendientes»): añade la **transcripción** del texto, el **nombre** del alumno
   (en «Datos privados») y **tapa la firma** si hace falta.

### Capturar notas nuevas
Al final de cada clase, abre el grupo del curso actual (fecha «Hoy» por defecto)
y fotografía las notas del día. Complétalas después igual que las históricas.

### Anexos (para no duplicar tu archivo docente)
En la pestaña **Anexos** de cada grupo puedes añadir:
- **Lista de clase** — valida los nombres y las fechas del capítulo. Siempre
  confidencial (solo versión con nombres); su foto no se muestra salvo que la
  reveles a propósito.
- **Actividades** — visibles en ambas versiones por defecto.
- **Valoraciones** — solo para el jefe por defecto.

### Componer y exportar el libro
1. En **⚙️ Ajustes** pon el título, subtítulo, autor, periodo y dedicatoria del
   libro, y elige una o dos notas por página.
2. En **📖 Tu libro** verás el libro montado (siempre anónimo).
3. «**Exportar a PDF**» abre el diálogo de impresión:
   - **Versión anónima** (recomendada): sin nombres, con las firmas tapadas.
   - **Versión con nombres — solo para el jefe**: marca la casilla de
     confirmación. Incluye nombres, fotos originales, lista de clase y
     valoraciones, con la marca «CONFIDENCIAL».
   En el destino de impresión elige **«Guardar como PDF»**.

---

## 4. Copias de seguridad

Tus datos viven **solo en este navegador**. Si borras los datos del navegador o
cambias de equipo, se pierden. Por eso:

- En **⚙️ Ajustes → Copias de seguridad**, usa «**Exportar todo**» con
  regularidad. Se descarga un archivo `.json` con todo (grupos, notas, anexos,
  ajustes y fotos). **Contiene nombres**: guárdalo en lugar seguro.
- «**Exportar solo ese grupo**» genera un archivo más pequeño, útil para pasar
  un capítulo del móvil al ordenador.
- «**Importar copia**» fusiona un archivo con lo que ya tengas (no duplica: cada
  elemento se actualiza por su identificador).
- Activa el **almacenamiento persistente** desde Ajustes para que el navegador
  no borre los datos si le falta espacio.

Los archivos de copia (`copia_libro_reconocimiento*.json`) y cualquier PDF que
guardes dentro de esta carpeta están excluidos del repositorio por `.gitignore`:
nunca se subirán por accidente.

---

## 5. Verificación

Prueba automatizada de humo (recorre el flujo completo y comprueba las garantías
de privacidad) con [Playwright](https://playwright.dev):

```bash
NODE_PATH="$(npm root -g)" node libro_reconocimiento/verificacion/prueba_humo.cjs
```

No forma parte de ninguna integración continua; es un script manual.

### Comprobación manual rápida
1. Abre `index.html`; con las herramientas de desarrollo en modo sin conexión,
   confirma que no hace ninguna petición de red.
2. Crea un grupo, importa una foto (mejor una foto **vertical** de móvil real:
   debe verse bien orientada).
3. Añade transcripción y un nombre; comprueba que el nombre no aparece en la
   galería ni en el visor.
4. Tapa una firma: la galería y el visor deben mostrar la versión tapada;
   reábrela y quita el rectángulo para comprobar que vuelve el original.
5. Exporta el PDF anónimo y el PDF con nombres desde Chrome («Guardar como
   PDF»): revisa que el anónimo no tiene nombres ni firmas y que no hay imágenes
   cortadas entre páginas.
6. Recarga el navegador: los datos siguen ahí y el libro está en anónimo.
