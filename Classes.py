import numpy as np
from matplotlib import pyplot as plt
PERIOD = 0.002

class Dot:
    def __init__(self, t, q, v):
        self.t = t
        self.q = q
        self.v = v

    def get_Dot_to_Dot_plot_line(Dot1, Dot2):
        if Dot1.t == Dot2.t:
            return [], [], [], []
        num_ticks = 1000
        t = np.linspace(Dot1.t, Dot2.t, num_ticks)
        a = (Dot2.v - Dot1.v)/(Dot2.t - Dot1.t)
        b = (Dot1.v*Dot2.t - Dot2.v*Dot1.t)/(Dot2.t - Dot1.t)

        q = Dot1.q + 0.5 * a * t**2 + b * t - \
            (0.5 * a * Dot1.t**2 + b * Dot1.t)
        v = a * t + b
        a = np.full(num_ticks, a)

        return t, q, v, a

    def print_Dots(Dots):
        template = "{0:>10}|{1:>10}|{2:>10}|{3:>10}"
        print("Dots:")
        print(template.format("№", "t", "q", "v"))
        template = "{0:10}|{1:10.3f}|{2:10.3f}|{3:10.3f}"
        for i, Dot in enumerate(Dots):
            print(template.format(i, Dot.t, Dot.q, Dot.v))


class Joint:
    def __init__(self, positions, dq_max, ddq_max, trajectory_junction_ratio=0, do_numerical_control=False):
        self.positions = positions
        self.dq_max = dq_max
        self.ddq_max = ddq_max
        self.trajectory_junction_ratio = trajectory_junction_ratio


        self.trajectories = []
        self.calculate_path(t0=0)
        if do_numerical_control:
            self.numerical_control()

    def calculate_trajectory(q0, qf, dq_m, ddq_m, t0):
        dq = abs(qf-q0)
        # triangular check
        c = np.sqrt(dq*ddq_m)

        if c <= dq_m:
            t1 = np.sqrt(dq/ddq_m)
            dq_m = ddq_m * t1
            T = t1
            tf = 2*t1
        else:
            t1 = dq_m/ddq_m
            T = dq/dq_m
            tf = T+t1

        dq_m *= np.sign(qf-q0)
        ddq_m *= np.sign(qf-q0)

        t0_Dot = Dot(t0, q0, 0)
        t1_Dot = Dot(t0 + t1, q0 + (0.5*ddq_m*t1**2),
                         dq_m)
        T_Dot = Dot(t0 + T, q0 + (0.5*ddq_m*t1**2) +
                        dq_m*(T-t1), dq_m)
        tf_Dot = Dot(t0 + tf, qf, 0)
        return t0_Dot, t1_Dot, T_Dot, tf_Dot

    def __recalculate_positions_after_junction(junction_Dot, trajectory, qf, dq_max, ddq_max):
        t0_Dot, t1_Dot, T_Dot, tf_Dot = trajectory
        t, q, v, a = Dot.get_Dot_to_Dot_plot_line(
            junction_Dot, t1_Dot)
        t1_Dot.q = q[-1]
        dq = qf - q[-1]
        dq_max = np.sign(dq) * abs(dq_max)
        ddq_max = np.sign(dq) * abs(ddq_max)

        T = t1_Dot.t + (dq - 0.5 * ddq_max *
                          (t1_Dot.t - t0_Dot.t)**2)/dq_max
        T_Dot.t = T
        T_Dot.q = q[-1] + (T - t1_Dot.t) * dq_max
        tf_Dot.t = T + t1_Dot.t - t0_Dot.t

    def calculate_path(self, t0=0):
        self.trajectories = []
        for i in range(1, len(self.positions)):
            t0_Dot, t1_Dot, T_Dot, tf_Dot = Joint.calculate_trajectory(
                self.positions[i-1], self.positions[i], self.dq_max, self.ddq_max, t0)

            # Here we will recalculate points 
            # for trajectory junction

            if self.trajectory_junction_ratio > 0 and len(self.trajectories):
                traj = [t0_Dot, t1_Dot, T_Dot, tf_Dot]
                Joint.__recalculate_positions_after_junction(
                    self.trajectories[-1][2], traj, self.positions[i], self.dq_max, self.ddq_max)

            self.trajectories.append([t0_Dot, t1_Dot, T_Dot, tf_Dot])
            t0 = T_Dot.t + (1-self.trajectory_junction_ratio) * \
                (tf_Dot.t - T_Dot.t)

    def numerical_control(self):
        # Getting points in correct order
        Dots = []
        Dots.append(self.trajectories[0][0])
        Dots.append(self.trajectories[0][1])
        Dots.append(self.trajectories[0][2])
        for i in range(1, len(self.trajectories)):
            Dots.append(self.trajectories[i][0])
            Dots.append(self.trajectories[i-1][3])
            Dots.append(self.trajectories[i][1])
            Dots.append(self.trajectories[i][2])
        Dots.append(self.trajectories[-1][3])

        # Recalculate times after numerical control
        prev_time, prev_time_modified = 0, 0
        for Dot in Dots:
            t_modified = prev_time_modified + \
                np.ceil((Dot.t - prev_time)/PERIOD) * PERIOD
            prev_time = Dot.t
            Dot.t = t_modified
            prev_time_modified = t_modified

        # Recalculate velocities after numerical control
        for i in range(len(self.trajectories)):
            trajectory = self.trajectories[i]
            new_v = (self.positions[i+1] - self.positions[i]
                     )/(trajectory[2].t - trajectory[0].t)
            trajectory[1].v = new_v
            trajectory[2].v = new_v

        # Recalculate positions after numerical control
        for i in range(len(self.trajectories)):
            trajectory = self.trajectories[i]
            dq_m = trajectory[1].v
            ddq_m = (trajectory[1].v - trajectory[0].v) / \
                (trajectory[1].t - trajectory[0].t)

            if self.trajectory_junction_ratio > 0 and i > 0:
                Joint.__recalculate_positions_after_junction(
                    self.trajectories[i-1][2], trajectory, self.positions[i+1], dq_m, ddq_m)
            else:
                trajectory[1].q = trajectory[0].q + \
                    (0.5*ddq_m*(trajectory[1].t - trajectory[0].t)**2)
                trajectory[2].q = trajectory[1].q + trajectory[1].v * \
                    (trajectory[2].t - trajectory[1].t)

    def get_path_Dots(self):
        Dots = [Dot for trajectory in self.trajectories for Dot in trajectory]
        if self.trajectory_junction_ratio == 0:
            return Dots

        path_Dots = [Dots[0]]
        for i in range(len(Dots)):
            if i % 4 == 1 or i % 4 == 2:
                path_Dots.append(Dots[i])
        path_Dots.append(Dots[-1])

        return path_Dots

    def __plot(self, X, Y, index, name, vlines=[], hlines=[]):
        plt.subplot(1, 3, index)
        plt.plot(X, Y, linewidth=2, label=name, c='r')
        plt.xlabel('t (s)', fontsize=18)
        plt.ylabel(f'{name}(t) ($\degree$/s)', fontsize=18)
        plt.grid(color='black', linestyle='--', linewidth=1.0, alpha=0.7)
        plt.grid(True)
        plt.xlim([min(X), max(X)])
        plt.ylim([min(Y), 1.1 * max(Y)])
        if vlines:
            plt.vlines(vlines, min(Y), 1.1 * max(Y),
                       linestyles='--', linewidth=2)
        if hlines:
            plt.hlines(hlines, min(Y), 1.1 * max(Y),
                       linestyles='--', linewidth=2)

    def plot_path_Dots(self, print_Dots=False):
        Dots = self.get_path_Dots()
        if print_Dots:
            Dot.print_Dots(Dots)
            # Dot.print_Dots([Dot for trajectory in self.trajectories for Dot in trajectory])

        T, Q, V, A = [], [], [], []
        for i in range(1, len(Dots)):
            t, q, v, a = Dot.get_Dot_to_Dot_plot_line(
                Dots[i-1], Dots[i])
            T = np.concatenate((T, t))
            Q = np.concatenate((Q, q))
            V = np.concatenate((V, v))
            A = np.concatenate((A, a))

        plt.figure(figsize=(18, 5))
        self.__plot(T, Q, 1, 'q')
        self.__plot(T, V, 2, 'v')
        self.__plot(T, A, 3, 'a')
        plt.show()
