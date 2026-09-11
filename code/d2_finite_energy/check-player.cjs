// Tests the animation's event handlers with a minimal DOM stub, not a browser.
const fs = require("node:fs");
const path = require("node:path");
const vm = require("node:vm");
const assert = require("node:assert/strict");
const html = fs.readFileSync(path.resolve(__dirname,
  "../../writing/d2_finite_energy/animations/core_collapse.html"), "utf8");
const scripts = [...html.matchAll(/<script\b([^>]*)>([\s\S]*?)<\/script>/g)];
const frames = JSON.parse(scripts.find(s => s[1].includes("frame-data"))[2]);
assert.equal(frames.length, 72);
assert.ok(frames.every(f => f.startsWith("data:image/webp;base64,")));
assert.equal(new Set(frames).size, 72);
const elements = {};
for (const id of ["frame", "position", "counter", "play", "reset", "speed", "frame-data"]) {
  elements[id] = {value: id === "speed" ? "1" : "0", handlers: {}, attrs: {},
    addEventListener(event, fn) { this.handlers[event] = fn; },
    setAttribute(key, val) { this.attrs[key] = val; }};
}
elements["frame-data"].textContent = JSON.stringify(frames);
const timers = new Map();
let next = 0;
const document = {hidden: false, handlers: {},
  getElementById(id) { return elements[id]; },
  addEventListener(event, fn) { this.handlers[event] = fn; }};
const context = vm.createContext({document,
  setInterval(fn) { const id = ++next; timers.set(id, fn); return id; },
  clearInterval(id) { timers.delete(id); }});
vm.runInContext(scripts.find(s => !s[1].includes("frame-data"))[2], context);
assert.equal(timers.size, 0); // no autoplay
assert.equal(elements.frame.src, frames[0]);
elements.play.handlers.click();
assert.equal(timers.size, 1);
[...timers.values()][0]();
assert.equal(elements.frame.src, frames[1]);
elements.play.handlers.click();
assert.equal(timers.size, 0);
elements.position.value = "50";
elements.position.handlers.input();
assert.equal(elements.frame.src, frames[50]);
elements.reset.handlers.click();
assert.equal(elements.frame.src, frames[0]);
elements.play.handlers.click();
elements.speed.value = "2";
elements.speed.handlers.change();
assert.equal(timers.size, 1);
for (let step = 0; step < 75 && timers.size; ++step) [...timers.values()][0]();
assert.equal(elements.frame.src, frames[71]);
assert.equal(timers.size, 0); // stops at endpoint
elements.play.handlers.click();
assert.equal(elements.frame.src, frames[0]); // replay
document.hidden = true;
document.handlers.visibilitychange();
assert.equal(timers.size, 0);
console.log("Animation event checks passed: no autoplay, play/pause, scrub, reset, speed, stop, replay, hidden pause.");
