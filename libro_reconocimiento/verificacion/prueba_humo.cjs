/*
 * Prueba de humo del Libro de reconocimiento (Playwright, sin dependencias del repo).
 *
 * Verifica el recorrido completo por file:// y, sobre todo, las garantías de
 * privacidad: que ningún nombre de alumno ni documento confidencial se filtra
 * a la versión anónima ni sobrevive en el DOM tras imprimir.
 *
 * Requisitos (ya presentes en el entorno de desarrollo):
 *   - Node 18+  ·  paquete «playwright» instalado (global u otro alcance)
 *   - Chromium disponible para Playwright
 *
 * Ejecución:
 *   NODE_PATH="$(npm root -g)" node libro_reconocimiento/verificacion/prueba_humo.cjs
 *
 * No forma parte de ningún CI: es un script manual de verificación.
 */
'use strict';
const path = require('path');
const os = require('os');
const fs = require('fs');
const zlib = require('zlib');

let chromium;
try { ({ chromium } = require('playwright')); }
catch (e) {
  console.error('No se encontró «playwright». Instálalo o define NODE_PATH="$(npm root -g)".');
  process.exit(2);
}

const RUTA_APP = 'file://' + path.resolve(__dirname, '..', 'index.html');

/* ---- PNG RGB de color sólido, generado sin dependencias ---- */
function crc32(buf) {
  let c = ~0;
  for (let i = 0; i < buf.length; i++) { c ^= buf[i]; for (let k = 0; k < 8; k++) c = (c >>> 1) ^ (0xEDB88320 & -(c & 1)); }
  return (~c) >>> 0;
}
function chunk(tipo, datos) {
  const t = Buffer.from(tipo, 'ascii');
  const len = Buffer.alloc(4); len.writeUInt32BE(datos.length, 0);
  const crc = Buffer.alloc(4); crc.writeUInt32BE(crc32(Buffer.concat([t, datos])), 0);
  return Buffer.concat([len, t, datos, crc]);
}
function makePng(w = 400, h = 300, rgb = [210, 150, 80]) {
  const firma = Buffer.from([137, 80, 78, 71, 13, 10, 26, 10]);
  const ihdr = Buffer.alloc(13);
  ihdr.writeUInt32BE(w, 0); ihdr.writeUInt32BE(h, 4);
  ihdr[8] = 8; ihdr[9] = 2;
  const fila = Buffer.alloc(1 + w * 3);
  for (let x = 0; x < w; x++) { fila[1 + x*3] = rgb[0]; fila[2 + x*3] = rgb[1]; fila[3 + x*3] = rgb[2]; }
  const idat = zlib.deflateSync(Buffer.concat(Array.from({ length: h }, () => fila)));
  return Buffer.concat([firma, chunk('IHDR', ihdr), chunk('IDAT', idat), chunk('IEND', Buffer.alloc(0))]);
}

const NOMBRE = 'NOMBRE_SECRETO_XYZ';
const LISTA = 'LISTA_SECRETA_XYZ';

(async () => {
  const navegador = await chromium.launch();
  const contexto = await navegador.newContext({ acceptDownloads: true });
  const page = await contexto.newPage();
  const errores = [];
  page.on('pageerror', e => errores.push('pageerror: ' + e.message));
  page.on('console', m => { if (m.type() === 'error') errores.push('console: ' + m.text()); });
  page.on('dialog', d => d.type() === 'prompt' ? d.accept('BORRAR') : d.accept());

  // Stub de impresión: captura el HTML impreso y dispara afterprint sin diálogo nativo.
  await page.addInitScript(() => {
    window.__impresiones = [];
    window.print = () => {
      window.__impresiones.push(document.getElementById('libro-impresion').innerHTML);
      window.dispatchEvent(new Event('afterprint'));
    };
  });

  let n = 0;
  const ok = (cond, msg) => { if (!cond) throw new Error('FALLO: ' + msg); n++; console.log('  ✔ ' + msg); };
  const buffer = makePng();

  console.log('Prueba de humo — Libro de reconocimiento\n' + RUTA_APP + '\n');
  await page.goto(RUTA_APP);
  await page.waitForSelector('#vista-grupos:not([hidden])');
  ok(await page.isVisible('[data-test="grupos-vacio"]'), 'arranca por file:// con estado vacío');

  // --- Grupo ---
  await page.click('[data-test="btn-nuevo-grupo"]');
  await page.fill('[data-test="input-grupo-nombre"]', 'Generación 2019');
  await page.fill('[data-test="input-grupo-curso"]', '2019-2020');
  await page.click('[data-test="btn-guardar-grupo"]');
  await page.click('[data-test="tarjeta-grupo"] .tarjeta-grupo-cuerpo');
  await page.waitForSelector('#vista-galeria:not([hidden])');

  // --- Captura por tanda con fecha aproximada (solo año) ---
  await page.selectOption('[data-test="captura-precision"]', 'anio');
  await page.fill('[data-test="captura-fecha-anio"]', '2019');
  await page.setInputFiles('[data-test="input-archivos"]', [
    { name: 'a.png', mimeType: 'image/png', buffer },
    { name: 'b.png', mimeType: 'image/png', buffer },
  ]);
  await page.waitForFunction(() => document.querySelectorAll('[data-test="nota-miniatura"]').length === 2);
  ok(true, 'captura por tanda: 2 fotos → 2 notas');
  ok((await page.textContent('#num-pendientes')) === '2', 'las notas de la tanda quedan pendientes');
  ok((await page.textContent('[data-test="nota-miniatura"] .sutil')).includes('2019'), 'fecha aproximada (2019) aplicada');

  // --- Rellenar una nota con nombre privado ---
  await page.click('[data-test="nota-miniatura"]:first-child');
  await page.waitForSelector('#vista-nota:not([hidden])');
  const notaId = await page.evaluate(() => location.hash.split('/')[2]);
  await page.fill('[data-test="nota-transcripcion"]', 'Gracias por todo, profe.');
  await page.click('[data-test="datos-privados"] summary');
  await page.fill('[data-test="nota-alumno"]', NOMBRE);

  // --- Tapar la firma en la foto ---
  await page.click('[data-test="btn-tapar"]');
  await page.waitForSelector('#vista-censura:not([hidden])');
  const caja = await page.locator('[data-test="lienzo-censura"]').boundingBox();
  await page.mouse.move(caja.x + caja.width * 0.2, caja.y + caja.height * 0.6);
  await page.mouse.down();
  await page.mouse.move(caja.x + caja.width * 0.8, caja.y + caja.height * 0.85, { steps: 8 });
  await page.mouse.up();
  await page.click('[data-test="btn-censura-guardar"]');
  await page.waitForSelector('#vista-nota:not([hidden])');
  const tieneAnon = await page.evaluate(id => new Promise(res => {
    const r = indexedDB.open('libro_reconocimiento');
    r.onsuccess = () => { const g = r.result.transaction('imagenes','readonly').objectStore('imagenes').get(id+':anonima'); g.onsuccess = () => res(!!g.result); };
  }), notaId);
  ok(tieneAnon, 'tapar nombres genera la variante anónima');
  await page.click('[data-test="btn-guardar-nota"]');
  await page.waitForSelector('#vista-galeria:not([hidden])');
  ok(!(await page.content()).includes(NOMBRE), 'el nombre del alumno NO aparece en la galería');

  // --- Anexo: lista de clase (solo jefe) ---
  await page.click('[data-test="pestana-anexos"]');
  await page.click('[data-test="btn-nuevo-anexo"]');
  await page.selectOption('[data-test="nuevo-anexo-tipo"]', 'lista_clase');
  await page.setInputFiles('[data-test="input-anexo-archivo"]', { name: 'lista.png', mimeType: 'image/png', buffer });
  await page.waitForSelector('#vista-anexo:not([hidden])');
  ok(await page.isHidden('[data-test="anexo-imagen"]'), 'la lista de clase no muestra su foto por defecto');
  await page.fill('[data-test="anexo-titulo"]', LISTA);
  await page.click('[data-test="btn-guardar-anexo"]');
  await page.waitForSelector('#vista-galeria:not([hidden])');

  // --- Libro (visor, siempre anónimo) ---
  await page.click('[data-test="btn-ir-libro"]');
  await page.waitForSelector('[data-test="libro-visor"] .pagina');
  const visor = await page.innerHTML('[data-test="libro-visor"]');
  ok(!visor.includes(NOMBRE), 'el visor NO muestra el nombre del alumno');
  ok(!visor.includes(LISTA), 'el visor NO muestra la lista de clase (solo jefe)');

  // --- Exportar anónimo ---
  await page.click('[data-test="btn-exportar-pdf"]');
  await page.waitForSelector('#dialogo-exportar[open]');
  ok(await page.isChecked('[data-test="export-anonimo"]'), 'el diálogo abre en anónimo por defecto');
  await page.click('[data-test="btn-abrir-impresion"]');
  await page.waitForFunction(() => window.__impresiones.length === 1);
  const anon = await page.evaluate(() => window.__impresiones[0]);
  ok(!anon.includes(NOMBRE) && !anon.includes(LISTA), 'PDF anónimo: sin nombre ni lista de clase');
  ok(anon.includes('pagina-portada'), 'PDF anónimo: con portada');

  // --- Exportar con nombres (requiere confirmación) ---
  await page.click('[data-test="btn-exportar-pdf"]');
  await page.waitForSelector('#dialogo-exportar[open]');
  await page.check('[data-test="export-nombres"]');
  await page.click('[data-test="btn-abrir-impresion"]');
  await page.waitForTimeout(250);
  ok((await page.evaluate(() => window.__impresiones.length)) === 1, 'sin confirmar el checkbox NO exporta con nombres');
  await page.check('[data-test="confirmar-nombres"]');
  await page.click('[data-test="btn-abrir-impresion"]');
  await page.waitForFunction(() => window.__impresiones.length === 2);
  const conNombres = await page.evaluate(() => window.__impresiones[1]);
  ok(conNombres.includes(NOMBRE), 'PDF con nombres: incluye el nombre del alumno');
  ok(conNombres.includes(LISTA), 'PDF con nombres: incluye la lista de clase (validación documental)');
  ok(conNombres.includes('CONFIDENCIAL'), 'PDF con nombres: lleva marca CONFIDENCIAL');
  await page.waitForFunction(() => !document.getElementById('libro-impresion').innerHTML.includes('NOMBRE_SECRETO_XYZ'));
  ok(true, 'tras imprimir, el DOM se limpia de nombres (afterprint)');

  // --- @media print oculta la app ---
  await page.emulateMedia({ media: 'print' });
  ok(!(await page.isVisible('#app')), 'en @media print, #app se oculta');
  await page.emulateMedia({ media: 'screen' });

  // --- Copia de seguridad: exportar → borrar → importar ---
  await page.click('[data-test="btn-ir-ajustes"]');
  await page.waitForSelector('#vista-ajustes:not([hidden])');
  const [descarga] = await Promise.all([
    page.waitForEvent('download'),
    page.click('[data-test="btn-exportar-todo"]'),
  ]);
  const ruta = path.join(os.tmpdir(), 'humo_backup_' + Date.now() + '.json');
  await descarga.saveAs(ruta);
  const backup = JSON.parse(fs.readFileSync(ruta, 'utf8'));
  ok(backup.formato === 'libro_reconocimiento_backup' && backup.esquema === 1, 'copia con formato y esquema válidos');
  ok(backup.notas[0].alumno_nombre === NOMBRE, 'la copia conserva el nombre (íntegra)');

  await page.click('[data-test="btn-borrar-todo"]');
  await page.waitForFunction(() => document.querySelectorAll('[data-test="tarjeta-grupo"]').length === 0);
  await page.click('[data-test="btn-ir-ajustes"]');
  await page.waitForSelector('#vista-ajustes:not([hidden])');
  await page.setInputFiles('[data-test="input-importar"]', ruta);
  await page.waitForFunction(() =>
    Array.from(document.querySelectorAll('.aviso-exito')).some(a => a.textContent.includes('Importado')), { timeout: 20000 });
  const restaurado = await page.evaluate(() => new Promise(res => {
    const r = indexedDB.open('libro_reconocimiento');
    r.onsuccess = () => { const c = r.result.transaction('notas','readonly').objectStore('notas').count(); c.onsuccess = () => res(c.result); };
  }));
  ok(restaurado === 2, 'importar restaura las notas');
  fs.unlinkSync(ruta);

  // --- Persistencia tras recarga: el libro vuelve a anónimo ---
  await page.goto(RUTA_APP + '#/libro');
  await page.waitForSelector('[data-test="libro-visor"] .pagina');
  const trasRecarga = await page.innerHTML('[data-test="libro-visor"]');
  ok(!trasRecarga.includes(NOMBRE) && !trasRecarga.includes('CONFIDENCIAL'), 'tras recargar, el libro está en anónimo');

  ok(errores.length === 0, 'sin errores de página/consola' + (errores.length ? ': ' + JSON.stringify(errores) : ''));

  await navegador.close();
  console.log(`\n✔ Prueba de humo superada (${n} comprobaciones).`);
})().catch(e => { console.error('\nERROR: ' + e.message); process.exit(1); });
