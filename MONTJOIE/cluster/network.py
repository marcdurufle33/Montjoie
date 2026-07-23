# Copyright (C) 2005-2006 Vivien Mallet
# 
# This files provides facilities to deal with computations over a Linux
# network.
# 
# This file is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
# 
# This file is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License (file ``license'') for more details.


"""Provides class Network designed to managed processes launched over the network."""

import os, sys, types, string, socket, subprocess, pwd

ssh_command = "ssh "
comp_host_tuples = [
    ("devel01", 24), \
    ("devel02", 24), \
    ("devel03", 24) ]


class Host:
    """Dedicated to host management."""
    def __init__(self, host = socket.gethostname()):
        if isinstance(host, tuple) or isinstance(host, list):
            self.name = str(host[0])
            if len(host) == 1:
                self.cpu = 1
            else:
                self.cpu = int(host[1])
        else:
            self.name = str(host)
            self.cpu = 1
    def LoadAverage(self, t = 60.):
        """Returns the load average within t seconds."""
        if self.name == socket.gethostname():
            return os.getloadavg()
        else:
            (s, o) = Network([self.name]).LaunchWait("__get_loc_loadavg", \
                                                     t, host = self.name)
            
            if s == 0:
                sortie = o.split()
                sortie = sortie[(len(sortie)-3):len(sortie)]		
                return tuple([float(elt) for elt in sortie ])
            else:
                return (9999, 9999, 9999)
    def LoadAverageFile(self):
        """Returns the load average within t seconds in a file."""
        import tempfile, warnings
        warnings.simplefilter("ignore")
        file = os.tempnam("/tmp/")
        if self.name == socket.gethostname():
            sfile = open(file, 'w')
            sfile.write(' '.join([str(x) for x in os.getloadavg()]))
            sfile.close()
            return [file, ""]
        else:
            sfile = open(file, 'w')
            sfile.close()
            cmd_file_name = os.tempnam("/tmp/")
            cmd_file = open(cmd_file_name, 'w')
            cmd_file.write("__get_loc_loadavg >& " + file)
            cmd_file.close()
            os.chmod(cmd_file.name, 0o700)
            Network([self.name]).LaunchScreen(cmd_file.name, host = self.name)
            return [file, cmd_file_name]
    def LaunchInt(self, command):
        """Launches a command in interactive mode (using os.system)."""
        return os.system(command)
    def LaunchFG(self, command):
        """Launches a command in the foreground."""
        return subprocess.getstatusoutput(command)
    def LaunchBG(self, command):
        """Launches a command in the background and returns a Popen4 object."""
        # Output is redirected.
        command = "( " + command + "; ) &> /dev/null"
        return subprocess.popen(command, shell=True, executable="/bin/bash")
    def LaunchWait(self, command, ltime, wait = 0.1):
        """Launches a command in the foreground and waits for its output""" \
        + """ for a given time after which the process is killed."""
        import time
        # Output is redirected.
        file_index = 0
        while os.path.exists("/tmp/output-" + str(file_index)):
            file_index += 1
        open("/tmp/output-" + str(file_index), 'w')
        file_name = "/tmp/output-" + str(file_index)
        os.chmod(file_name, 0o600)
        command = "( " + command + "; ) &> " + file_name
        p = subprocess.popen(command, shell=True, executable="/bin/bash")
        t = time.time()
        while (p.poll() == None and time.time() - t < ltime):
            time.sleep(wait)
        s = p.poll()
        if s == None:
            try:
                os.kill(p.pid, 9)
            except:
                pass
        o = open("/tmp/output-" + str(file_index), 'r').read()
        os.remove(file_name)
        return (s, o)
        

class Network:
    """Dedicated to network management."""
    def __init__(self, hosts = []):
        if hosts == []:
            self.hosts = []
            return
        if isinstance(hosts, str):
            hosts = comp_host_tuples
            
        self.hosts = []
        for host in hosts:
            self.hosts.append(Host(host))
    def PrintHostNames(self):
        for host in self.hosts: print(host.name)
    def GetHostNames(self):
        l = []
        for host in self.hosts: l.append(host.name)
        return l
    def LaunchInt(self, command, host = socket.gethostname()):
        """Launches a command in interactive mode (using os.system)."""
        if not isinstance(host, Host):
            host = Host(host)
        if host.name == socket.gethostname():
            return os.system(command)
        else:
            return os.system(ssh_command + host.name + " \"" + command + "\"")
    def LaunchFG(self, command, host = socket.gethostname()):
        """Launches a command in the foreground."""
        if not isinstance(host, Host):
            host = Host(host)
        if host.name == socket.gethostname():
            (s, o) = subprocess.getstatusoutput(command)
            # # "commands.getstatusoutput" removes the last '\n' if it exists!
            # # Assuming that there was a line break:
            # o += "\n"
        else:
            (s, o) = subprocess.getstatusoutput(ssh_command + host.name + " \"" + command + "\"")
        # Note: "commands.getstatusoutput" removes the last '\n' if it exists!
        return (s, o)
    def LaunchWait(self, command, ltime, wait = 0.1, host = socket.gethostname()):
        """Launches a command in the foreground and waits for its output""" \
        + """ for a given time after which the process is killed."""
        import time
        if not isinstance(host, Host):
            host = Host(host)
        # Output is redirected.
        file_index = 0
        while os.path.exists("/tmp/output-" + str(file_index)):
            file_index += 1
        open("/tmp/output-" + str(file_index), 'w')
        file_name = "/tmp/output-" + str(file_index)
        os.chmod(file_name, 0o600)
        if host.name == socket.gethostname():
            command = "( " + command + "; ) &> " + file_name
        else:
            command = "( "+ssh_command + host.name + " \"" + command + "\"; ) &> " + file_name
        p = subprocess.Popen(command,shell=True, executable="/bin/bash")
        t = time.time()
        while (p.poll() == None and time.time() - t < ltime):
            time.sleep(wait)
        s = p.poll()
        if s == None:
            try:
                os.kill(p.pid, 9)
            except:
                pass
        o = open("/tmp/output-" + str(file_index), 'r').read()
        os.remove(file_name)
        return (s, o)
    def LaunchBG(self, command, host = socket.gethostname()):
        """Launches a command in the background and returns a Popen4 object."""
        if not isinstance(host, Host):
            host = Host(host)
        # Output is redirected.
        command = "( " + command + "; ) "
        if host.name == socket.gethostname():
            return subprocess.Popen(command,shell=True, executable="/bin/bash")
        else:
            return subprocess.Popen(ssh_command + host.name + " \"" + command + "\"", shell=True, executable="/bin/bash")
    def LaunchScreen(self, command, host = socket.gethostname()):
        """Launches a command (file) in a screen."""
        if not isinstance(host, Host):
            host = Host(host)
        if host.name != socket.gethostname():
            file = open(command)
            lines = file.readlines()
            file.close()
            file = open(command, 'w')
            file.write(ssh_command + host.name + " " +'\n'.join(lines))
            file.close()
        command = "screen -d -m \"" + command + "\""
        os.system(command)
    def SendMail(self, subject, msg = "", toaddr = pwd.getpwuid(os.getuid())[0], \
                 fromaddr = pwd.getpwuid(os.getuid())[0]):
        from email.MIMEText import MIMEText
        if toaddr.find("@") == -1:
            toaddr += "@cerea.enpc.fr"
        if fromaddr.find("@") == -1:
            fromaddr += "@cerea.enpc.fr"
        msg = MIMEText(msg)
        msg['Subject'] = subject
        msg['From'] = fromaddr
        msg['To'] = toaddr
        import smtplib
        server = smtplib.SMTP('localhost')
        server.sendmail(fromaddr, toaddr, msg.as_string())
        server.quit()
    def SendMailAttach(self, subject, attachments, msg = "",
                       toaddr = pwd.getpwuid(os.getuid())[0], \
                       fromaddr = pwd.getpwuid(os.getuid())[0]):
        from email.MIMEMultipart import MIMEMultipart
        from email.MIMEText import MIMEText
        if toaddr.find("@") == -1:
            toaddr += "@cerea.enpc.fr"
        if fromaddr.find("@") == -1:
            fromaddr += "@cerea.enpc.fr"
        tmp = MIMEMultipart()
        tmp.attach(MIMEText(msg))
        msg = tmp
        msg['Subject'] = subject
        msg['From'] = fromaddr
        msg['To'] = toaddr
        if type(attachments) == types.StringType:
            file = open(attachments, 'rb')
            msg.attach(MIMEText(file.read()))
            file.close()
        else:
            for attachment in attachments:
                file = open(attachment, 'rb')
                msg.attach(MIMEText(file.read()))
                file.close()
        import smtplib
        server = smtplib.SMTP('localhost')
        server.sendmail(fromaddr, toaddr, msg.as_string())
        server.quit()
    def GetLoadAverages(self, time = 5.):
        l = []
        for host in self.hosts:
            l.append(host.LoadAverage(time))
        return l
    def GetAvailableHost(self, load_limit = 0.3, time = 5.):
        import time
        load_min = 100000
        host_min = self.hosts[0]
        files = []
        for host in self.hosts:
            files.append(host.LoadAverageFile())
        time.sleep(2)
        for ihost in range(len(self.hosts)):
            file = open(files[ihost][0])
            lines = file.readlines()
            file.close()
            if len(lines) != 0:
                avr = lines[0].split()[0]
                try:
                    avr = float(avr)
                except:
                    avr = 9999
            else:
                avr = 9999
            l = avr
            if l < load_limit:
                for jhost in range(len(self.hosts)):
                    for file in files[jhost]:
                        if file != "":
                            os.remove(file)
                return self.hosts[ihost].name
            else:
                if l < load_min:
                    load_min = l
                    host_min = self.hosts[ihost]
        for ihost in range(len(self.hosts)):
            for file in files[ihost]:
                if file != "":
                    os.remove(file)
        return host_min.name
    def GetRandomAvailableHost(self, load_limit = 0.3):
        host_available = []
        for host in self.hosts:
            l = host.LoadAverage()[0]
            if l < load_limit:
                host_available.append(host)
        if host_available:
            import random
            return host_available[random.randint(0, len(host_available) - 1)].name
        else:
            return self.hosts[0].name

if __name__ == "__main__":
    f = Network("comp")
    # f.SendMailAttach("Essai", ["test", "test"], "La prochaine fois, j'attache un fichier..")
    print(f.GetAvailableHost())
    print(f.GetHostNames())
