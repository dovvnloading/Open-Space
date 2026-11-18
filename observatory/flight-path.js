// Draw a great-circle slew on an equirectangular equatorial sky chart.
// The small difference between ICRS and J2000 axes is below the displayed precision.
const rad = Math.PI / 180;
function vector(c) {
  const a = c.raDeg * rad,
    b = c.decDeg * rad;
  return [Math.cos(b) * Math.cos(a), Math.cos(b) * Math.sin(a), Math.sin(b)];
}
export function pointingArc(from, to) {
  const a = vector(from),
    b = vector(to);
  const angle = Math.acos(
    Math.min(
      1,
      Math.max(
        -1,
        a.reduce((s, x, i) => s + x * b[i], 0),
      ),
    ),
  );
  const degrees = angle / rad;
  if (angle < 1e-7 || Math.PI - angle < 1e-7) return { degrees, path: "" };
  const segments = Math.max(2, Math.ceil(degrees / 2));
  const path = [];
  let lastX;
  for (let i = 0; i <= segments; i++) {
    const u = Math.sin((1 - i / segments) * angle) / Math.sin(angle);
    const v = Math.sin((i / segments) * angle) / Math.sin(angle);
    const [x, y, z] = a.map((value, k) => u * value + v * b[k]);
    const sx = 360 - ((Math.atan2(y, x) / rad + 360) % 360);
    const sy = 90 - Math.asin(Math.min(1, Math.max(-1, z))) / rad;
    const begin = lastX === undefined || Math.abs(lastX - sx) > 180;
    path.push((begin ? "M" : "L") + sx.toFixed(2) + " " + sy.toFixed(2));
    lastX = sx;
  }
  return { degrees, path: path.join(" ") };
}
