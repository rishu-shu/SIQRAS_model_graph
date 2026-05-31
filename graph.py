import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider, Button

# =========================
# Time settings
# =========================
dt = 0.1
T = 50
time = np.arange(0, T, dt)

# =========================
# Initial conditions
# =========================
S0, I0, Q0, R0, A0 = 90, 5, 0, 0, 5

# =========================
# Default parameters
# =========================
beta0 = 0.01      # infection rate
gamma0 = 0.05     # quarantine rate
alphaSA0 = 0.01   # S → A
alphaIA0 = 0.02   # I → A

# Fixed parameters
delta = 0.1       # Q → R
epsilon = 0.03    # I → R
sigma = 0.01      # loss of immunity
phi = 0.02        # quarantine failure

# =========================
# SIQRAS model (N=0, mu=0)
# =========================
def simulate(beta, gamma, alphaSA, alphaIA):
    S = np.zeros(len(time))
    I = np.zeros(len(time))
    Q = np.zeros(len(time))
    R = np.zeros(len(time))
    A = np.zeros(len(time))

    S[0], I[0], Q[0], R[0], A[0] = S0, I0, Q0, R0, A0

    for t in range(1, len(time)):
        dS = - beta*S[t-1]*I[t-1] - alphaSA*S[t-1]*A[t-1] + sigma*R[t-1]
        dI = beta*S[t-1]*I[t-1] - alphaIA*I[t-1]*A[t-1] - (epsilon + gamma)*I[t-1]
        dQ = gamma*I[t-1] - (phi + delta)*Q[t-1]
        dR = delta*Q[t-1] - sigma*R[t-1] + epsilon*I[t-1]
        dA = alphaIA*I[t-1]*A[t-1] + alphaSA*S[t-1]*A[t-1]

        S[t] = S[t-1] + dS*dt
        I[t] = I[t-1] + dI*dt
        Q[t] = Q[t-1] + dQ*dt
        R[t] = R[t-1] + dR*dt
        A[t] = A[t-1] + dA*dt

    return S, I, Q, R, A

# =========================
# Initial simulation
# =========================
S, I, Q, R, A = simulate(beta0, gamma0, alphaSA0, alphaIA0)

# =========================
# Plot
# =========================
fig, ax = plt.subplots()
plt.subplots_adjust(bottom=0.35)

lS, = ax.plot(time, S, label="Susceptible (S)")
lI, = ax.plot(time, I, label="Infected (I)")
lQ, = ax.plot(time, Q, label="Quarantined (Q)")
lR, = ax.plot(time, R, label="Recovered (R)")
lA, = ax.plot(time, A, label="Antidotal (A)")

ax.set_title("SIQRAS Malware Spread Model (N=0,μ=0)")
ax.set_xlabel("Time")
ax.set_ylabel("Number of Systems")
ax.legend()
ax.grid(True)

# =========================
# Sliders
# =========================
ax_beta = plt.axes([0.15, 0.25, 0.7, 0.03])
ax_gamma = plt.axes([0.15, 0.20, 0.7, 0.03])
ax_alphaSA = plt.axes([0.15, 0.15, 0.7, 0.03])
ax_alphaIA = plt.axes([0.15, 0.10, 0.7, 0.03])

s_beta = Slider(ax_beta, 'β Infection', 0.0, 0.03, valinit=beta0)
s_gamma = Slider(ax_gamma, ' γ Quarantine', 0.0, 0.6, valinit=gamma0)
s_alphaSA = Slider(ax_alphaSA, 'αSA (S→A)', 0.0, 0.05, valinit=alphaSA0)
s_alphaIA = Slider(ax_alphaIA, 'αIA (I→A)', 0.0, 0.05, valinit=alphaIA0)

def update(val):
    S, I, Q, R, A = simulate(
        s_beta.val,
        s_gamma.val,
        s_alphaSA.val,
        s_alphaIA.val
    )
    lS.set_ydata(S)
    lI.set_ydata(I)
    lQ.set_ydata(Q)
    lR.set_ydata(R)
    lA.set_ydata(A)
    fig.canvas.draw_idle()

s_beta.on_changed(update)
s_gamma.on_changed(update)
s_alphaSA.on_changed(update)
s_alphaIA.on_changed(update)

# =========================
# WannaCry preset
# =========================
ax_wc = plt.axes([0.35, 0.02, 0.3, 0.05])
btn_wc = Button(ax_wc, "WannaCry Preset")

def wannacry(event):
    s_beta.set_val(0.025)     # very high infection
    s_gamma.set_val(0.04)     # delayed quarantine
    s_alphaSA.set_val(0.005)  # weak early patching
    s_alphaIA.set_val(0.01)   # slow antivirus cleanup

btn_wc.on_clicked(wannacry)

plt.show()
