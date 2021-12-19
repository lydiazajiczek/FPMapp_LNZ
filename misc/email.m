function email(subject,message)
% Summary of this function goes here
%   Detailed explanation goes here
props = java.lang.System.getProperties;
props.setProperty('mail.smtp.port','587');
props.setProperty('mail.smtp.auth','true');
props.setProperty('mail.smtp.starttls.enable','true');
setpref('Internet','E_mail','christopher.bendkowski.18@ucl.ac.uk');
setpref('Internet','SMTP_Username', 'apikey');
setpref('Internet','SMTP_Password', getenv('sendgridkey'));
setpref('Internet','SMTP_Server','smtp.sendgrid.net');
sendmail('christopher.bendkowski.18@ucl.ac.uk',subject,message);
end

