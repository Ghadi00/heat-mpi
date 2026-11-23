# Network setup, 2026-11-23

All VMs: VirtualBox, NAT networking, Tailscale inside the VM.

## Cross-VM ping-pong runs 
- heat-lynn: coffee-shop public Wi-Fi (shared with heat-rim)
- heat-rim:  same coffee-shop public Wi-Fi
- heat-ghadi: home Wi-Fi (different network)

Note: lynn and rim were on the same Wi-Fi, but each VM sits behind VirtualBox NAT
and public Wi-Fi usually isolates clients, so their traffic went out through the
cafe's internet connection and back. Both directions shared one throttled uplink,
which matches rim<->lynn having the highest latency and lowest bandwidth.

## Same-VM runs, path check (paths.txt)
All three VMs on different networks (each at their own location).
- paths.txt therefore describes this later setup, not the cross-VM runs.
