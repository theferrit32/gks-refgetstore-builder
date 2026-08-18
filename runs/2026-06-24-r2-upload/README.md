# 2026-06-24 R2 upload

This is a **reconstructed** record assembled from the former `UPLOAD_PLAN.md`
and `TODO.md`. The raw upload log is unavailable, so the recorded outcome is
limited to facts preserved in those notes.

- Destination: `theferrit32-public:theferrit32-public/refgetstore/2026-06-24`
- Command: `./upload_store.sh --bucket theferrit32-public --transfers 64 --checkers 64`
- Concurrency: 64 transfers and 64 checkers (16 was abandoned as too slow on
  the tiny-object tail)
- Uploaded objects: 783,641
- Verification: `rclone check` reported 783,641 matching files and 0 differences
- Public URL: <https://static.ferriter.dev/refgetstore/2026-06-24/>
- Viewer: <https://refget.databio.org/explore-store/overview?url=https%3A%2F%2Fstatic.ferriter.dev%2Frefgetstore%2F2026-06-24%2F>

The notes also recorded that the scoped R2 token could read and write objects
but could not list buckets. No credentials are retained here.
