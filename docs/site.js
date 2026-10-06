'use strict';
const copy = document.querySelector('[data-copy-note]');
if (copy) copy.addEventListener('click', async () => {
  const note = document.getElementById('submission-note').textContent;
  try { await navigator.clipboard.writeText(note); copy.textContent = 'Copied'; }
  catch (_) { copy.textContent = 'Select and copy the note below'; }
});
const feed = document.getElementById('feed-status');
if (feed) fetch('data/source_health.json').then(r => { if (!r.ok) throw new Error('Feed unavailable'); return r.json(); }).then(data => {
  const item = data.sources.find(s => s.id === 'leaderboard');
  if (item && item.status === 'parsed' && typeof item.top_score === 'number') {
    feed.textContent = `Latest automated leaderboard check: ${item.top_score.toFixed(4)} · ${item.checked_at.slice(0, 16).replace('T',' ')} UTC. Reachability is not validation of our model.`;
  } else {
    feed.textContent = `Automated check ${data.checked_at.slice(0,10)} could not refresh the leaderboard. The 2026-10-06 research snapshot remains 0.3345; do not treat it as live.`;
  }
}).catch(() => { feed.textContent = 'Live feed unavailable. Displaying the dated 2026-10-06 research snapshot, not a live score.'; });
